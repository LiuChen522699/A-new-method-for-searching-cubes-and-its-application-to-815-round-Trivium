#include "Trivium.h"


// --------------------------------------------- 利用代数次数估计筛项 ----------------------------------------------
Poly& FilterTerm(char Infile[], Ivterm &CubeFlag, int16_t OutPutRound, uint8_t IvCon[], char Outfile[])
{
    // 从文件中读入多项式
    const int dim = CubeFlag.count();
    int16_t TargetRound = 0;
    static Poly ExpandingANF(10);
    ReadANF(TargetRound, ExpandingANF, Infile);

    // 获得对应cube在指定轮数下各个状态的代数次数
    int State[288] = { 0 };
    int SState[288] = { 0 };

    int16_t TempRound = OutPutRound - TargetRound;
    cout << "Round = " << TempRound << endl;
    Numeric_Mapping(State, SState, CubeFlag, TempRound);

    // 利用给定的代数次数进行筛项
    Poly* TempANF = new Poly(10);
    TempANF->SetPolyLen(ExpandingANF.Size);
    for (uint32_t pt = 0; pt < ExpandingANF.Size; pt++)
    {
        int TestDeg = 0;
        for (uint16_t pd = 0; pd < 288; pd++)
            if ((ExpandingANF.poly[pt].pterm[(pd >> 5)] >> (31 - (pd & 0x1f))) & 1)
            {
                if ((pd != 287) && ((ExpandingANF.poly[pt].pterm[((pd + 1) >> 5)] >> (31 - ((pd + 1) & 0x1f))) & 1))
                    TestDeg += SState[pd++];
                else
                    TestDeg += State[pd];
            }
        if (TestDeg >= dim)
            TempANF->poly[TempANF->Size++] = ExpandingANF.poly[pt];
    }
    Poly TmpANF(*TempANF);
    delete TempANF;
    cout << "Save Terms(Pass Numeric Mapping): " << TmpANF.Size << endl;

    ExpandingANF.SetPolyLen(TmpANF.Size);

clock_t starttime = clock();
// -------------------------------------------------------------------------------------------------
#if 1  // 直接利用MILP筛项    (300-Round: 422.185s)
    for (uint32_t pd = 0; pd < TmpANF.Size; pd++)
    {
        if (TriviumEval(CubeFlag, TmpANF.poly[pd].pterm, IvCon, TempRound))
            ExpandingANF.poly[ExpandingANF.Size++] = TmpANF.poly[pd];
        cout << "Round: " << TempRound << ", All Terms: " << TmpANF.Size << ", Have Tested Terms: " << pd+1 << ", Accepted Terms: " << ExpandingANF.Size << endl;
    }
#endif

#if 0  // 边利用MILP筛项，边利用整除丢弃无用项    (300-Round:  444.471s)
    Poly* DiscardTerms = new Poly(10);
    DiscardTerms->SetPolyLen(TmpANF.Size);
    ExpandingANF.SetPolyLen(TmpANF.Size);
    uint32_t pd = 0;
    while (pd < TmpANF.Size)
    {
        int degflag = TmpANF.poly[pd].deg;
        // 每次对最大次数的项进行估计，并对利用无用项对后面低次的项进行筛选，看是否进行估计
        while ((TmpANF.poly[pd].deg == degflag) && (pd < TmpANF.Size))
        {
            if (TmpANF.poly[pd].deg != 0)
            {
                if (TriviumEval(CubeFlag, TmpANF.poly[pd].pterm, IvCon, TempRound))
                    ExpandingANF.poly[ExpandingANF.Size++] = TmpANF.poly[pd];
                else
                    DiscardTerms->poly[DiscardTerms->Size++] = TmpANF.poly[pd];
                cout << "All Terms: " << TmpANF.Size << ", Have Tested Terms: " << pd+1 
                     << ", Accepted Terms: " << ExpandingANF.Size << ", Discard Terms: " << DiscardTerms->Size << endl;
            }
            pd += 1;
        }

        uint32_t UselessNum = 0;
        if (DiscardTerms->Size > 0)
            for (uint32_t pt = pd; pt < TmpANF.Size; pt++)
                for (uint32_t pf = 0; pf < DiscardTerms->Size; pf++)
                    if (Divisibility(DiscardTerms->poly[pf], TmpANF.poly[pt]))
                    {
                        UselessNum += 1;
                        TmpANF.poly[pt].deg = 0;
                        break;
                    }
        cout << "UselessNum = " << UselessNum << endl;
        DiscardTerms->Size = 0;
    }
    delete DiscardTerms;
#endif 
// -------------------------------------------------------------------------------------------------
clock_t endtime = clock();

    cout << "Save Terms(Pass CBDP_MILP): " << ExpandingANF.Size << ", time = " << double(endtime - starttime) / CLOCKS_PER_SEC << endl;
    ExpandingANF.write_value(TempRound, Outfile, "Terms: " + to_string(ExpandingANF.Size));
    return ExpandingANF;
}

// InPutRound:表示要开始回推的轮数
Poly& FilterPartTerm(int16_t& InPutRound, char file[], Ivterm &CubeFlag, uint8_t IvCon[], uint32_t ConstrTermNum)
{
    // 从文件中读入多项式
    const int dim = CubeFlag.count();
    Poly* InPutANF  = new Poly;
    Poly* OutPutANF = new Poly;
    Poly UpdateFunc[8], TempPoly1(SIZE), TempPoly2(SIZE);
    Poly SNull(0), S1(4), S94(4), S178(4), S1S94(16), S1S178(16), S94S178(16), SAll(64);
    InvExpressInit(TempPoly1, S1, S94, S178, S1S94, S1S178, S94S178, SAll);
    UpdateFunc[0] = SNull;
    UpdateFunc[1] = S178;
    UpdateFunc[2] = S94;
    UpdateFunc[3] = S94S178;
    UpdateFunc[4] = S1;
    UpdateFunc[5] = S1S178; 
    UpdateFunc[6] = S1S94;
    UpdateFunc[7] = SAll;
    ReadANF(InPutRound, SNull, file);
    // InPutANF->SetPolyLen(100 * SIZE);
    // OutPutANF->SetPolyLen(100 * SIZE);
    InPutANF->PolyCopy(SNull);
    SNull.SetPolyLen(0);
    uint32_t SaveTermNum = 1;
    
    // SaveTermNum = 0, 表示已经没有剩余的项了, 因此不再需要回推了; InPutRound = 0, 表示已经回推到第0轮了, 因此不再需要回推
    while ( (SaveTermNum <= 2000) && (SaveTermNum) && (InPutRound) )
    {
        // ConstrTermNum: 表示一次检测的数量，最好不超过10000项
        while ((OutPutANF->Size <= 3000) && (InPutRound))
        {
            TempPoly1.Size = 0;
            TempPoly2.Size = 0;

            // 将输入多项式区分为新回推的部分和平移的部分
            for (uint32_t pt = 0; pt < InPutANF->Size; pt++)
            {
                uint8_t flag = (((InPutANF->poly[pt].pterm[0] >> 31) & 1) << 2) | (((InPutANF->poly[pt].pterm[2] >> 2) & 1) << 1) | ((InPutANF->poly[pt].pterm[5] >> 14) & 1);
                uint8_t flagweight = ((InPutANF->poly[pt].pterm[0] >> 31) & 1) + ((InPutANF->poly[pt].pterm[2] >> 2) & 1) + ((InPutANF->poly[pt].pterm[5] >> 14) & 1);
                for (int i = 0; i < 8; i++)
                    InPutANF->poly[pt].pterm[i] = (InPutANF->poly[pt].pterm[i] << 1) | ((InPutANF->poly[pt].pterm[i + 1] >> 31) & 1);
                InPutANF->poly[pt].pterm[8] <<= 1;
                InPutANF->poly[pt].deg -= flagweight;
                if (flag)
                {
                    InPutANF->poly[pt].pterm[2] &= 0xfffffff7;
                    InPutANF->poly[pt].pterm[5] &= 0xffff7fff;
                    PolyMul(TempPoly2, UpdateFunc[flag], InPutANF->poly[pt]);
                }
                else
                    TempPoly1.poly[TempPoly1.Size++] = InPutANF->poly[pt];
            }
            // 对新回推的项再向前展开或者平移
            for (uint32_t pt = 0; pt < OutPutANF->Size; pt++)
            {
                uint8_t flag = (((OutPutANF->poly[pt].pterm[0] >> 31) & 1) << 2) | (((OutPutANF->poly[pt].pterm[2] >> 2) & 1) << 1) | ((OutPutANF->poly[pt].pterm[5] >> 14) & 1);
                uint8_t flagweight = ((OutPutANF->poly[pt].pterm[0] >> 31) & 1) + ((OutPutANF->poly[pt].pterm[2] >> 2) & 1) + ((OutPutANF->poly[pt].pterm[5] >> 14) & 1);
                for (int i = 0; i < 8; i++)
                    OutPutANF->poly[pt].pterm[i] = (OutPutANF->poly[pt].pterm[i] << 1) | ((OutPutANF->poly[pt].pterm[i + 1] >> 31) & 1);
                OutPutANF->poly[pt].pterm[8] <<= 1;
                OutPutANF->poly[pt].deg -= flagweight;
                if (flag)
                {
                    OutPutANF->poly[pt].pterm[2] &= 0xfffffff7;
                    OutPutANF->poly[pt].pterm[5] &= 0xffff7fff;
                    PolyMul(TempPoly2, UpdateFunc[flag], OutPutANF->poly[pt]);
                }
                else
                    TempPoly2.poly[TempPoly2.Size++] = OutPutANF->poly[pt];
            }

            InPutANF->PolyCopy(TempPoly1);   // 里边是平移过来不需要回推的项，因此对于检测来说没有影响
            OutPutANF->PolyCopy(TempPoly2);  // 里边是要进行检测的项，项数达到一定数量就不再回推了
            OutPutANF->RemoveDup();
            InPutRound -= 1;
            cout << "OutPutANF Terms: " << (OutPutANF->Size + InPutANF->Size) << ", InPutRound = " << InPutRound << endl;
        }

        // 获得对应cube在指定轮数下各个状态的代数次数
        int State[288] = { 0 };
        int SState[288] = { 0 };
        Numeric_Mapping(State, SState, CubeFlag, InPutRound);

        // 利用给定的代数次数进行筛项
        TempPoly1.Size = 0;
        for (uint32_t pt = 0; pt < OutPutANF->Size; pt++)
        {
            int TestDeg = 0;
            for (uint16_t pd = 0; pd < 288; pd++)
                if ((OutPutANF->poly[pt].pterm[(pd >> 5)] >> (31 - (pd & 0x1f))) & 1)
                {
                    if ((pd != 287) && ((OutPutANF->poly[pt].pterm[((pd + 1) >> 5)] >> (31 - ((pd + 1) & 0x1f))) & 1))
                        TestDeg += SState[pd++];
                    else
                        TestDeg += State[pd];
                }
            if (TestDeg >= dim)
                TempPoly1.poly[TempPoly1.Size++] = OutPutANF->poly[pt];
        }
        cout << "Save Terms(Pass Numeric Mapping): " << (TempPoly1.Size + InPutANF->Size) << endl;

        OutPutANF->Size = 0;
        for (uint32_t pd = 0; pd < TempPoly1.Size; pd++)
        {
            if (TriviumEval(CubeFlag, TempPoly1.poly[pd].pterm, IvCon, InPutRound))
                InPutANF->poly[InPutANF->Size++] = TempPoly1.poly[pd];
            cout << "Round: " << InPutRound << ", All Terms: " << TempPoly1.Size << ", Have Tested Terms: " << pd << ", Accepted Terms: " << InPutANF->Size << endl;
        }
        InPutANF->RemoveDup();

        cout << "Save Terms(Pass CBDP_MILP): " << InPutANF->Size << endl;
        InPutANF->write_value(InPutRound, file, "Terms: " + to_string(InPutANF->Size));
        SaveTermNum = InPutANF->Size;
    }

    static Poly ExpandingANF;
    ExpandingANF = *InPutANF;
    delete InPutANF;
    delete OutPutANF;
    return ExpandingANF;
}


// ----------------------------------------- 针对给定项估计代数次(MILP CBDP) -----------------------------------------
bool TriviumEval(Ivterm& Cube, uint32_t Term[], uint8_t fflag[], int Round)
{
    
    GRBEnv *env = 0;
    env = new GRBEnv();
    GRBModel Model = GRBModel(*env);
    string Model_name = "Cube_Attack_to_Trivium_" + to_string(Round);  // 生成模型的名称
    Model.set(GRB_StringAttr_ModelName, Model_name.c_str());  // 模型命名
    Model.set(GRB_IntParam_LogToConsole, 0);
    Model.set(GRB_IntParam_Threads, 6);

    // Trivium 的更新的三个bit所需的抽头
    int Loc1[5] = { 65, 170, 90, 91, 92 };
    int Loc2[5] = { 161, 263, 174, 175, 176 };
    int Loc3[5] = { 242, 68, 285, 286, 287 };
    
    uint8_t iv_flag[80] = { 0 };  // 生成IV的flag
    uint8_t flag[80] = { 0 };  // 生成Cube的标志位,并补全iv_flag
    for (int i = 0; i < 80; i++)
    {
        if (Cube.test(i))
        {
            flag[i] = 1;
            iv_flag[i] = 2;
        }
        else
            iv_flag[i] = (fflag[i / 8] >> (7 - (i % 8))) & 0x1;
    }
    
    uint8_t VarFlag[288] = {0};  // 定义初态中IV对应的Flag
    string s_init_str[288];  // 添加初态, 首先定义生成变量所需的参数, 然后生成变量
    char s_type[288];
    double LB[288] = {0};
    double UB[288];
    for (int i = 0; i < 288; i++)
    {
        s_init_str[i] = "s_init_" + to_string(i);
        s_type[i] = 'B';
        UB[i] = 1;
    }
    GRBVar* s_state = Model.addVars(LB, UB, LB, s_type, s_init_str, 288);

    // 将最大代数次数最为求解目标函数
    GRBLinExpr key_sum;
    for (int i = 0; i < 80; i++)
        key_sum += s_state[i];
    Model.setObjective(key_sum, GRB_MAXIMIZE);

    // 将Cube变元添加到里边去
    for (int i = 0; i < 80; i++)
    {
        VarFlag[i] = 2;
        VarFlag[93 + i] = iv_flag[i];
        if (flag[i])
            Model.addConstr(s_state[93 + i] == 1, "Iv_constr_" + to_string(i));
        else
            Model.addConstr(s_state[93 + i] == 0, "Iv_constr_" + to_string(i));
    }
    for (int i = 80; i < 93; i++)
        Model.addConstr(s_state[i] == 0, "constant_constr_" + to_string(i));
    for (int i = 173; i < 288; i++)
        Model.addConstr(s_state[i] == 0, "constant_constr_" + to_string(i));
    VarFlag[285] = VarFlag[286] = VarFlag[287] = 1;

    // 更新迭代
    for (int r = 0; r < Round; r++)
    {
        TriviumCore(Model, s_state, VarFlag, Loc1, r);
        TriviumCore(Model, s_state, VarFlag, Loc2, r);
        TriviumCore(Model, s_state, VarFlag, Loc3, r);
        GRBVar temp = s_state[287];
        uint8_t tempflag = VarFlag[287];
        for (int i = 287; i > 0; i--)
        {
            s_state[i] = s_state[i-1];
            VarFlag[i] = VarFlag[i-1];
        }
        s_state[0] = temp;
        VarFlag[0] = tempflag;
    }

    // 增加一个条件, 给定一个项, 求该项的代数次数是否为dim
    GRBLinExpr TermConstr;
    uint32_t ZeroFlag = 0;
    for (int pd = 0; pd < 288; pd++)
    {
        if ((Term[pd >> 5] >> (31 - (pd & 0x1f))) & 1)
        {
            TermConstr += s_state[pd];
            if (VarFlag[pd] == 0)
                ZeroFlag = 1;
        }
        else
            Model.addConstr(s_state[pd] == 0, "FinalRoundState" + to_string(pd) + " to 0");
    }
    Model.addConstr(TermConstr >= 1, "TermIsExist");
    if (ZeroFlag)
        Model.addConstr(TermConstr == 0, "TermIsNotExist");
    
    // 模型求解, 并返回是否有解
    Model.optimize();
    if (Model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        delete env;
        return true;
    }
    else
    {
        delete env;
        return false;
    }
}

void TriviumCore(GRBModel &Model, GRBVar* Var, uint8_t VarFlag[], int loc[], int round)
{
    uint8_t temp_Flag[10];
    string var_name[10];
    char var_type[10] = {'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B'};
    double LB[10] = { 0 };
    double UB[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    for (int i = 0; i < 5; i++)
        var_name[i] = "y_" + to_string(round) + "_" + to_string(loc[i]);
    for (int i = 5; i < 9; i++)
        var_name[i] = "z_" + to_string(round) + "_" + to_string(i - 4);
    var_name[9] = "a_" + to_string(round);
    GRBVar* var_new = Model.addVars(LB, UB, LB, var_type, var_name, 10);

    for (int i = 0; i < 4; i++)
    {
        Model.addConstr(var_new[i] == Var[loc[i]] - var_new[5 + i], "Constr_" + to_string(round) + "_" + to_string(i));
        temp_Flag[i] = temp_Flag[5 + i] = VarFlag[loc[i]];
    }

    Model.addConstr(var_new[9] >= var_new[7], "Constr_" + to_string(round) + "_4");
    Model.addConstr(var_new[9] >= var_new[8], "Constr_" + to_string(round) + "_5");
    temp_Flag[9] = temp_Flag[7] * temp_Flag[8];

    if (temp_Flag[9] == 0)
        Model.addConstr(var_new[9] == 0, "Constr_" + to_string(round) + "_7");
    else if (temp_Flag[9] > 2)
        temp_Flag[9] = 2;

    Model.addConstr(var_new[4] == Var[loc[4]] + var_new[9] + var_new[5] + var_new[6], "Constr_" + to_string(round) + "_6");
    if ((VarFlag[loc[4]] > 1) || (temp_Flag[9] > 1) || (temp_Flag[5] > 1) || (temp_Flag[6] > 1))
        temp_Flag[4] = 2;
    else
        temp_Flag[4] = VarFlag[loc[4]] ^ temp_Flag[9] ^ temp_Flag[5] ^ temp_Flag[6];

    for (int i = 0; i < 5; i++)
    {
        Var[loc[i]] = var_new[i];
        VarFlag[loc[i]] = temp_Flag[i];
    }
}


// ----------------------------------------- 针对给定项估计代数次(数值映射) -----------------------------------------
// cube变元要求{0~79}, 给出指定轮数状态bit的代数次数
void Numeric_Mapping(int State[], int SState[], Ivterm& Cube, int Round)
{
    int dA[1152+93] = {0};
    int dB[1152+84] = {0};
    int dC[1152+111] = {0};

    int dAA[111] = { 0 };
    int dBB[111] = { 0 };
    int dCC[111] = { 0 };

    // 初始化
    for (int i = 80; i < 93; i++)
        dA[92 - i] = -9999;
    for (int i = 0; i < 84; i++)
        dB[83 - i] = -9999;
    for (int i = 0; i < 108; i++)
        dC[110 - i] = -9999;
    // 添加cube变元
    for (int i = 0; i < 80; i++)
        if (Cube.test(i))
            dB[83 - i] = 1;

    for (int t = 1; t <= Round; t++)
    {
        int lA[3] = {dA[(t-1)+27], dA[t-1], dB[(t-1)+6]};
        int lB[3] = {dB[(t-1)+15], dB[t-1], dC[(t-1)+24]};
        int lC[3] = {dC[(t-1)+45], dC[t-1], dA[(t-1)+24]};

        int dA_arr[2] = {Deg_Mul(dA, dB, dC, t, 0), max(lC, 3)};
        dA[93 + t - 1] = max(dA_arr, 2);
        int dB_arr[2] = {Deg_Mul(dA, dB, dC, t, 1), max(lA, 3)};
        dB[84 + t - 1] = max(dB_arr, 2);
        int dC_arr[2] = {Deg_Mul(dA, dB, dC, t, 2), max(lB, 3)};
        dC[111 + t - 1] = max(dC_arr, 2);

        if (t == Round)
        {
            dAA[0] = dB_arr[0];
            dBB[0] = dC_arr[0];
            dCC[0] = dA_arr[0];
        }
    }
    for (int i = 0; i < 93; i++)
        State[i] = dA[Round + 92 - i];
    for (int i = 0; i < 84; i++)
        State[93 + i] = dB[Round + 83 - i];
    for (int i = 0; i < 111; i++)
        State[177 + i] = dC[Round + 110 - i];

    for (int t = Round+1; t <= Round+110; t++)
    {
        int lA[3] = { dA[(t - 1) + 27], dA[t - 1], dB[(t - 1) + 6] };
        int lB[3] = { dB[(t - 1) + 15], dB[t - 1], dC[(t - 1) + 24] };
        int lC[3] = { dC[(t - 1) + 45], dC[t - 1], dA[(t - 1) + 24] };

        int dA_arr[2] = { Deg_Mul(dA, dB, dC, t, 0), max(lC, 3) };
        dA[93 + t - 1] = max(dA_arr, 2);
        int dB_arr[2] = { Deg_Mul(dA, dB, dC, t, 1), max(lA, 3) };
        dB[84 + t - 1] = max(dB_arr, 2);
        int dC_arr[2] = { Deg_Mul(dA, dB, dC, t, 2), max(lB, 3) };
        dC[111 + t - 1] = max(dC_arr, 2);

        dAA[t - Round] = dB_arr[0];
        dBB[t - Round] = dC_arr[0];
        dCC[t - Round] = dA_arr[0];
    }

    for (int i = 0; i < 92; i++)
        SState[91 - i] = dAA[i];
    for (int i = 0; i < 83; i++)
        SState[93 + 82 - i] = dBB[i];
    for (int i = 0; i < 110; i++)
        SState[177 + 109 - i] = dCC[i];
    SState[92] = State[92] + State[93];
    SState[176] = State[176] + State[177];
}


// 利用flag表示使用的是 (Deg_Mul_A, Deg_Mul_B, Deg_Mul_C) = (0, 1, 2)
int Deg_Mul(int dA[], int dB[], int dC[], int t, int flag)
{
    int *dX1, *dX2;
    if (flag == 0)
    {
        dX1 = dC;
        dX2 = dB;
    }
    else if (flag == 1)
    {
        dX1 = dA;
        dX2 = dC;
    }
    else if (flag == 2)
    {
        dX1 = dB;
        dX2 = dA;
    }
    else 
        return 0;

    int nn[3] = {93, 84, 111};
    int pho[3] = {2, 0, 1};

    int pho1 = pho[flag];
    int pho2 = pho[pho1];
    int t1 = t - nn[pho1] + 1;
    
    if (t1 <= 0)
        return dX1[nn[pho1] + t1 - 1] + dX1[nn[pho1] + t1];
    
    int t2 = t1 - nn[pho2] + 2;
    int temp_arr[3] = {dX2[nn[pho2] + t2 - 1] + dX1[nn[pho1] + t1],
                       dX2[nn[pho2] + t2 + 1] + dX1[nn[pho1] + t1 - 1],
                       dX2[nn[pho2] + t2 - 1] + dX2[nn[pho2] + t2] + dX2[nn[pho2] + t2 + 1]};

    int linear_f0[3][3] = {{dC[t1+45], dC[t1], dA[t1+24]}, 
                           {dA[t1+27], dA[t1], dB[t1+6]},
                           {dB[t1+15], dB[t1], dC[t1+24]}};
    int linear_f1[3][3] = {{dC[(t1-1)+45], dC[t1-1], dA[(t1-1)+24]},
                           {dA[(t1-1)+27], dA[t1-1], dB[(t1-1)+6]},
                           {dB[(t1-1)+15], dB[t1-1], dC[(t1-1)+24]}};

    int max_arr[3] = {min(temp_arr, 3), 
                      max(linear_f0[pho1], 3) + dX1[nn[pho1] + t1 - 1], 
                      max(linear_f1[pho1], 3) + dX1[nn[pho1] + t1]};
                      
    int d = max(max_arr, 3);

    return d;
}

int min(int arr[], int len)
{
    if (len == 0)
        return 0;
    int temp = arr[0];
    for (int i = 1; i < len; i++)
        if (temp > arr[i])
            temp = arr[i];
    return temp;
}

int max(int arr[], int len)
{
    if (len == 0)
        return 0;
    int temp = arr[0];
    for (int i = 1; i < len; i++)
        if (temp < arr[i])
            temp = arr[i];
    return temp;
}


// -------------------------------------------------- 多项式的回推 --------------------------------------------------

// 回推并返回ANF: OutPutRound: 表示输出的轮数; TargetRound: 表示要回推的轮数, 即将输出函数回推TargetRound轮, 但表达式中的状态是(OutPutRound - TargetRound)的状态
void ExpressRecursivelyForTrivium(const int16_t OutPutRound, const int16_t TargetRound, char file[])
{
    // 定义回推的更新函数
    Poly UpdateFunction[8], ZR(6);
    Poly SNull(0), S1(4), S94(4), S178(4), S1S94(16), S1S178(16), S94S178(16), SAll(64);
    InvExpressInit(ZR, S1, S94, S178, S1S94, S1S178, S94S178, SAll);
    UpdateFunction[0] = SNull;
    UpdateFunction[1] = S178;
    UpdateFunction[2] = S94;
    UpdateFunction[3] = S94S178;
    UpdateFunction[4] = S1;
    UpdateFunction[5] = S1S178; 
    UpdateFunction[6] = S1S94;
    UpdateFunction[7] = SAll;

    Poly* OutPutANF = new Poly;
    Poly* InPutANF  = new Poly;
    Poly* TempPolyPt;
    Poly TempPoly(100*SIZE);
    (*OutPutANF).SetPolyLen(100*SIZE);
    (*InPutANF).SetPolyLen(100*SIZE);
    (*InPutANF).PolyCopy(ZR);

    int16_t CountRound = OutPutRound;  // OutPutRound的输出实际上使用的是 (OutPutRound-1) 轮的状态
    int16_t FinalRound = OutPutRound - TargetRound + 1;  // 441, 841-401+1

    while((CountRound--) > FinalRound)
    {
        clock_t start = clock();
        ExpressOneRound(*OutPutANF, *InPutANF, UpdateFunction, TempPoly);
        clock_t end = clock();
        TempPolyPt = InPutANF;
        InPutANF = OutPutANF;
        OutPutANF = TempPolyPt;
        cout << "CountRound = " << CountRound << ", time = " << double(end - start) / CLOCKS_PER_SEC << endl;
    }
    ExpressOneRound(*OutPutANF, *InPutANF, UpdateFunction, TempPoly);
    // OutPutANF->write_value(TargetRound, file, "OutPutANF");
    OutPutANF->write_output(TargetRound, file, "OutPutANF");
    delete OutPutANF;
    delete InPutANF;
}

// 多项式的回推: 逐轮回推
void ExpressOneRound(Poly& OutPutANF, Poly& InPutANF, Poly UpdateFunc[], Poly& TempPoly)
{
    OutPutANF.Size = 0;
    for (uint32_t pt = 0; pt < InPutANF.Size; pt++)
    {
        uint8_t flag = ( ((InPutANF.poly[pt].pterm[0]>>31)&1) << 2) | ( ((InPutANF.poly[pt].pterm[2]>>2)&1) << 1 ) | ( (InPutANF.poly[pt].pterm[5]>>14)&1 );
        uint8_t flagweight = ((InPutANF.poly[pt].pterm[0]>>31)&1) + ((InPutANF.poly[pt].pterm[2]>>2)&1) + ( (InPutANF.poly[pt].pterm[5]>>14)&1 );
        for (int i = 0; i < 8; i++)
            InPutANF.poly[pt].pterm[i] = (InPutANF.poly[pt].pterm[i] << 1) | ( (InPutANF.poly[pt].pterm[i+1] >> 31) & 1);
        InPutANF.poly[pt].pterm[8] <<= 1;
        InPutANF.poly[pt].deg -= flagweight;

        if (flag)
        {
            InPutANF.poly[pt].pterm[2] &= 0xfffffff7;
            InPutANF.poly[pt].pterm[5] &= 0xffff7fff;
            PolyMul(OutPutANF, UpdateFunc[flag], InPutANF.poly[pt]);
        }
        else
            OutPutANF.poly[OutPutANF.Size++] = InPutANF.poly[pt];
    }

    OutPutANF.RemoveDup();
}


// 初始化输出函数, 以及更新函数, 并给出更新函数的乘积
void InvExpressInit(Poly& ZR, Poly& S1, Poly& S94, Poly& S178, Poly& S1S94, Poly& S1S178, Poly& S94S178, Poly& SAll)
{
    Poly TempMul(128);
    S1.poly[0].deg = S94.poly[0].deg = S178.poly[0].deg = 2;
    S1.poly[1].deg = S94.poly[1].deg = S178.poly[1].deg = 1;
    S1.poly[2].deg = S94.poly[2].deg = S178.poly[2].deg = 1;
    S1.poly[3].deg = S94.poly[3].deg = S178.poly[3].deg = 1;
    S1.Size = S94.Size = S178.Size = 4;

    S1.poly[0].pterm[8] = (1 << (31 ^ 285&0x1f)) ^ (1 << (31 ^ 286&0x1f));    // 286 287 //0x6;
    S1.poly[1].pterm[8] = (1 << (31 ^ 287&0x1f));                             // 288
    S1.poly[2].pterm[7] = (1 << (31 ^ 242&0x1f));                             // 243
    S1.poly[3].pterm[2] = (1 << (31 ^ 68 &0x1f));                              // 69

    S94.poly[0].pterm[2] = (1 << (31 ^ 90 &0x1f)) ^ (1 << (31 ^ 91&0x1f));     // 91 92
    S94.poly[1].pterm[2] = (1 << (31 ^ 92 &0x1f));                             // 93
    S94.poly[2].pterm[2] = (1 << (31 ^ 65 &0x1f));                             // 66
    S94.poly[3].pterm[5] = (1 << (31 ^ 170&0x1f));                            // 171

    S178.poly[0].pterm[5] = (1 << (31 ^ 174&0x1f)) ^ (1 << (31 ^ 175&0x1f));  // 175 176
    S178.poly[1].pterm[5] = (1 << (31 ^ 176&0x1f));                           // 177
    S178.poly[2].pterm[5] = (1 << (31 ^ 161&0x1f));                           // 162
    S178.poly[3].pterm[8] = (1 << (31 ^ 263&0x1f));                           // 264

    // 计算乘积
    PolyMul(S1S94, TempMul, S1, S94);
    PolyMul(S1S178, TempMul, S1, S178);
    PolyMul(S94S178, TempMul, S94, S178);
    PolyMul(SAll, TempMul, S1, S94, S178);

    ZR.poly[0].pterm[2] = (1 << (31 ^ 65 & 0x1f));   // 66  (2 1 1073741824)
    ZR.poly[1].pterm[2] = (1 << (31 ^ 92 & 0x1f));   // 93  (2 28 8)
    ZR.poly[2].pterm[5] = (1 << (31 ^ 161 & 0x1f));  // 162 (5 1 1073741824)
    ZR.poly[3].pterm[5] = (1 << (31 ^ 176 & 0x1f));  // 171 (5 16 32768)
    ZR.poly[4].pterm[7] = (1 << (31 ^ 242 & 0x1f));  // 243 (7 18 8192)
    ZR.poly[5].pterm[8] = (1 << (31 ^ 287 & 0x1f));  // 288 (8 31 1)
    ZR.poly[0].deg = ZR.poly[1].deg = ZR.poly[2].deg = ZR.poly[3].deg = ZR.poly[4].deg = ZR.poly[5].deg = 1;
    ZR.Size = 6;

#if 0  // 检查初始化的结果
    char file[] = "Init.txt";
    ZR.write_output(0, file, "ZR");
    S1.write_output(0, file, "S1");
    S94.write_output(0, file, "S94");
    S178.write_output(0, file, "S178");
    S1S94.write_output(0, file, "S1S94");
    S1S178.write_output(0, file, "S1S178");
    S94S178.write_output(0, file, "S94S178");
    SAll.write_output(0, file, "SAll");
#endif
}

void ReadANF(int16_t &TargetRound, Poly& InPutANF, char file[])
{
    fstream fin;
    string line1, c1;
    uint32_t temp = 0;

    fin.open(file, ios_base::in);
    getline(fin, line1);
    fin >> c1 >> TargetRound >> temp;
    cout << line1 << endl << c1 << ' ' << TargetRound << ' ' << temp << endl;
    InPutANF.SetPolyLen(temp);
    InPutANF.Size = temp;

    clock_t starttime = clock();
    for (uint32_t pt = 0; pt < InPutANF.Size; pt++)
    {
        for (uint8_t pd = 0; pd < 9; pd++)
            fin >> hex >> InPutANF.poly[pt].pterm[pd];
        fin >> hex >> temp;
        InPutANF.poly[pt].deg = temp;
    }
        
    clock_t endtime = clock();
    cout << "Load ANF: " << double(endtime - starttime) / CLOCKS_PER_SEC << " s" << endl;

    fin.close();
    return ;
}


// --------------------------------------------- 对剩余的项利用三子集恢复 ----------------------------------------------

PolyKey& recover_supperpoly_BDPT(Poly& InPutANF, Ivterm& Cube, int Round, uint8_t iv_con[], char file[])
{
    AllVars Cube_term(0), Iv_term(0);
    for (int i = 0; i < 80; i++)
    {
        if ((iv_con[i >> 3] >> (7 ^ 7&i)) & 1)
            Iv_term.set(80 + i);
        if (Cube.test(i))
            Cube_term.set(80 + i);
    }
    Iv_term |= Cube_term;
    Iv_term.flip();

    PolyKey Superpoly;
    for (uint32_t pt = 0; pt < InPutANF.Size; pt++)
        bool flag = TriviumEvalBDPT(Superpoly, InPutANF.poly[pt].pterm, Cube_term, Iv_term, Round);
    sort(Superpoly.begin(), Superpoly.end(), CMP);
    // 接下里对获得的Superpoly进行排序, 然后整理超多项式, 并给出结果
    static PolyKey RetPoly;
    uint32_t pt = 0;
    for (; pt < Superpoly.size()-1; pt++)
        if (Superpoly[pt] == Superpoly[pt+1])
            pt++;
        else
            RetPoly.push_back(Superpoly[pt]);
    
    write2file(RetPoly, file, Cube_term);

    return RetPoly;
}

bool TriviumEvalBDPT(PolyKey& ResPoly, uint32_t Pterm[], AllVars &Cube, AllVars &Iv_term, int Round)
{
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_PoolSearchMode, 2);  // do a systematic search for the k- best solutions
    env.set(GRB_IntParam_PoolSolutions, 2000000000);  // Limit how many solutions to collect
    env.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);  // Limit the search space by setting a gap for the worst possible solution that will be acce
    env.set(GRB_IntParam_Threads, 10);
    GRBModel Model = GRBModel(env);  // 创建模型
    Model.set(GRB_StringAttr_ModelName, "Cube_Attack_to_trivium_" + to_string(Round) + " Round");
    // Model.set(GRB_DoubleParam_TimeLimit, 3600 * 24);  // 每次最多跑4h

    // Trivium 的更新的三个bit所需的抽头
    int Loc1[5] = { 65, 170, 90, 91, 92 };
    int Loc2[5] = { 161, 263, 174, 175, 176 };
    int Loc3[5] = { 242, 68, 285, 286, 287 };

    // 添加初态, 首先定义生成变量所需的参数, 然后生成变量
    string s_init_str[288];
    char s_type[288];
    double LB[288] = {0};
    double UB[288];
    for (int i = 0; i < 288; i++)
    {
        s_init_str[i] = "s_init_" + to_string(i);
        s_type[i] = 'B';
        UB[i] = 1;
    }
    GRBVar* s_state = Model.addVars(LB, UB, LB, s_type, s_init_str, 288);

    GRBLinExpr key_sum;
    for (int i = 0; i < 80; i++)
        key_sum += s_state[i];
    Model.setObjective(key_sum, GRB_MAXIMIZE);

    // 将Cube变元添加到里边去
    for (int i = 0; i < 80; i++)
    {
        if (Cube.test(80 + i))
            Model.addConstr(s_state[93 + i] == 1, "Iv_constr_" + to_string(i));
        else if (Iv_term.test(80 + i))
            Model.addConstr(s_state[93 + i] == 0, "Iv_constr_" + to_string(i));
    }
    for (int i = 80; i < 93; i++)
        Model.addConstr(s_state[i] == 0, "constant_constr_" + to_string(i));
    for (int i = 173; i < 285; i++)
        Model.addConstr(s_state[i] == 0, "constant_constr_" + to_string(i));

    // 更新迭代
    for (int r = 0; r < Round; r++)
    {
        TriviumCoreBDPT(Model, s_state, Loc1, r);
        TriviumCoreBDPT(Model, s_state, Loc2, r);
        TriviumCoreBDPT(Model, s_state, Loc3, r);
        GRBVar temp = s_state[287];
        for (int i = 287; i > 0; i--)
            s_state[i] = s_state[i - 1];
        s_state[0] = temp;
    }

    // 缺一个限制条件, 直接恢复一个给定的项
    GRBLinExpr TermConstr;
    for (int pd = 0; pd < 288; pd++)
    {
        if ((Pterm[pd >> 5] >> (31 - (pd & 0x1f))) & 1)
            TermConstr += s_state[pd];
        else
            Model.addConstr(s_state[pd] == 0, "FinalRoundState" + to_string(pd) + " to 0");
    }
    Model.addConstr(TermConstr >= 1, "TermIsExist");

    // 求解模型
    Model.optimize();
    
    if (Model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        int nSolutions = Model.get(GRB_IntAttr_SolCount);
        cout << " Number of solutions found : " << nSolutions << endl;
        PolyKey TempPoly, TmpPoly;
        TermCounter termcounter;

        for (int res = 0; res < nSolutions; res++)
        {
            AllVars p_term(0);
            Model.set(GRB_IntParam_SolutionNumber, res);

            for (int loc = 0; loc < 80; loc++)
            {
                GRBVar TempVar1 = Model.getVarByName("s_init_" + to_string(loc));
                GRBVar TempVar2 = Model.getVarByName("s_init_" + to_string(loc + 93));

                if (TempVar1.get(GRB_DoubleAttr_Xn) >= 0.5)
                    p_term.set(loc);
                if (TempVar2.get(GRB_DoubleAttr_Xn) >= 0.5)
                    p_term.set(80 + loc);
            }
            p_term ^= Cube;
            TempPoly.push_back(p_term);
        }

        // 处理超多项式的重复问题
        if (TempPoly.size() > 1)
            sort(TempPoly.begin(), TempPoly.end(), CMP);

        TmpPoly.push_back(TempPoly[0]);
        termcounter.push_back(1);
        for (int loc = 1; loc < TempPoly.size(); loc++)
        {
            if (*(TmpPoly.end() - 1) == TempPoly[loc])
                *(termcounter.end() - 1) += 1;
            else
            {
                TmpPoly.push_back(TempPoly[loc]);
                termcounter.push_back(1);
            }
        }
        TempPoly.clear();

        for (int loc = 1; loc < TmpPoly.size(); loc++)
            if (termcounter[loc]&1)
                ResPoly.push_back(TmpPoly[loc]);
        return true;
    }
    else 
        return false;
}

void TriviumCoreBDPT(GRBModel &Model, GRBVar* Var, int loc[], int round)
{
    GRBVar* var_new;
    string var_name[10];
    char var_type[10] = { 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B' };
    double UB[10] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    double LB[10] = { 0 };

    for (int i = 0; i < 5; i++)
        var_name[i] = "y_" + to_string(round) + '_' + to_string(loc[i]);
    for (int i = 5; i < 9; i++)
        var_name[i] = "z_" + to_string(round) + '_' + to_string(loc[i]) + '_' + to_string(i - 4);
    var_name[9] = "a_" + to_string(round) + '_' + to_string(loc[0]);
    var_new = Model.addVars(LB, UB, LB, var_type, var_name, 10);

    // Model for Copy
    for (int i = 0; i < 4; i++)
    {
        Model.addConstr(Var[loc[i]] <= var_new[i] + var_new[5 + i], "Copy1 Constr_" + to_string(round) + '_' + to_string(i));
        Model.addConstr(Var[loc[i]] >= var_new[5 + i], "Copy2 Constr_" + to_string(round) + '_' + to_string(i));
        Model.addConstr(Var[loc[i]] >= var_new[i], "Copy3 Constr_" + to_string(round) + '_' + to_string(i));
    }

    // Model for And
    Model.addConstr(var_new[9] == var_new[7], "And1 Constr_" + to_string(round) + "_4");
    Model.addConstr(var_new[9] == var_new[8], "And2 Constr_" + to_string(round) + "_4");

    // Model for Xor
    Model.addConstr(var_new[4] == Var[loc[4]] + var_new[9] + var_new[5] + var_new[6], "Xor Constr_" + to_string(round) + "_5");

    // replace variables 
    for (int i = 0; i < 5; i++)
        Var[loc[i]] = var_new[i];
}

bool CMP(const AllVars &p1, const AllVars &p2)
{
    if (p1.count() > p2.count())
        return true;
    else if (p1.count() < p2.count())
        return false;

    for (int i = 0; i < 160; i++)
        if (p1[i] > p2[i])
            return true;
        else if (p1[i] < p2[i])
            return false;
    return true;
}

void write2file(PolyKey &Superpoly, char file[], AllVars Cube_term)
{
    fstream fout;
    // 写入所求得的模型求得的结果
    fout.open(file, ios_base::app);
    fout << "Cube: ";
    for (int i = 0; i < 80; i++)
        if (Cube_term.test(80 + i))
            fout << 'v' << i << ' ';
    fout << endl << "Superpoly = ";

    int term_flag = 0;
    if (Superpoly.size())
    {
        for (int loc = 0; loc < Superpoly.size(); loc++)
        {
            if (term_flag)
                fout << " + ";
            term_flag += 1;

            if (Superpoly[loc] == 0)
                fout << '1';
            else
                for (int i = 0; i < 160; i++)
                    if (Superpoly[loc][i] && (i < 80))
                        fout << 'x' << i;
                    else if (Superpoly[loc][i] && (i >= 80))
                        fout << 'v' << i - 80;
        }
        if (!term_flag)
            fout << 0;
        fout << endl << "********************************************************" << endl;
        fout.close();
    }
    else
        fout << 0 << endl << "Superpoly Size: " << term_flag << "\n********************************************************\n";
}


// ---------------------------------------- 正向计算ANF, 将160个Key和Iv作为变元 -----------------------------------------
void Trivium(PolyFw S1[], PolyFw S2[], PolyFw S3[], Ivterm &CubeFlag, int Round, char file_c[])
{
    string files = file_c;
    files = files.substr(0, files.length() - 4) + '(' + to_string(Round) + ").txt";
    char file[20];
    files.copy(file, files.length(), 0);
    *(file + files.length()) = '\0';

    for (int i = 0; i < 80; i++)
    {
        TermFw temp1, temp2;
        temp1.pterm[(80 + i) >> 5] = (1 << (0x1f ^ ((80 + i) & 0x1f)));
        temp1.deg = 0;
        S1[i].poly[0] = temp1;
        S1[i].Size = 1;

        if (CubeFlag.test(i))
        {
            temp2.pterm[i >> 5] = (1 << (0x1f ^ (i & 0x1f)));
            temp2.deg = 1;
            S2[i].poly[0] = temp2;
            S2[i].Size = 1;
        }
    }

    TermFw temp;
    S3[110].poly[0] = S3[109].poly[0] = S3[108].poly[0] = temp;
    S3[110].Size = S3[109].Size = S3[108].Size = 1;

#if 0 // 输出初态
    string message;
    for (uint32_t pt = 0; pt < 93; pt++)
    {
        message = "  " + to_string(pt);
        S1[pt].write_poly(0, file, message);
    }
    for (uint32_t pt = 0; pt < 84; pt++)
    {
        message = "  " + to_string(93 + pt);
        S2[pt].write_poly(0, file, message);
    }
    for (uint32_t pt = 0; pt < 111; pt++)
    {
        message = "  " + to_string(177 + pt);
        S3[pt].write_poly(0, file, message);
    }
#endif

#if 1
    uint8_t BlockSize = 32;
    uint8_t SubCount = Round / BlockSize;
    uint8_t RestRound = Round - SubCount * BlockSize;
    cout << "RestRound = " << int(RestRound) << endl;

    PolyFw* temps1 = new PolyFw[BlockSize];
    PolyFw* temps2 = new PolyFw[BlockSize];
    PolyFw* temps3 = new PolyFw[BlockSize];

    PolyFw* bigtemp  = new PolyFw(SIZEFw * 20);
    PolyFw* temp_mul = new PolyFw(SIZEFw * 20);
    PolyFw* t1temp = new PolyFw;
    PolyFw* t2temp = new PolyFw;

    for (uint8_t count = 0; count < SubCount; count++)
    {
        for (uint8_t kk = 0; kk < BlockSize; kk++)
        {
            PolyMul(*bigtemp, *temp_mul, S1[90 - kk], S1[91 - kk]);
            PolyAdd(*t1temp, S1[65 - kk], S1[92 - kk]);
            PolyAdd(*t2temp, S2[77 - kk], *t1temp);
            PolyAdd(temps1[kk], *bigtemp, *t2temp);

            PolyMul(*bigtemp, *temp_mul, S2[81 - kk], S2[82 - kk]);
            PolyAdd(*t1temp, S2[68 - kk], S2[83 - kk]);
            PolyAdd(*t2temp, S3[86 - kk], *t1temp);
            PolyAdd(temps2[kk], *bigtemp, *t2temp);

            PolyMul(*bigtemp, *temp_mul, S3[108 - kk], S3[109 - kk]);
            PolyAdd(*t1temp, S3[kk + 45], S3[110 - kk]);
            PolyAdd(*t2temp, S1[68 - kk], *t1temp);
            PolyAdd(temps3[kk], *bigtemp, *t2temp);
        }

        for (int pd = 92; pd >= BlockSize; pd--)
            S1[pd].PolyCopy(S1[pd - BlockSize]);
        for (int pd = 83; pd >= BlockSize; pd--)
            S2[pd].PolyCopy(S2[pd - BlockSize]);
        for (int pd = 110; pd >= BlockSize; pd--)
            S3[pd].PolyCopy(S3[pd - BlockSize]);

        for (int pd = 0; pd < BlockSize; pd++)
        {
            S1[pd].PolyCopy(temps3[BlockSize - 1 - pd]);
            S2[pd].PolyCopy(temps1[BlockSize - 1 - pd]);
            S3[pd].PolyCopy(temps2[BlockSize - 1 - pd]);
        }
        cout << "Round : " << (count+1) * BlockSize << endl;
    }

    if (RestRound > 0)
    {
        for (uint8_t kk = 0; kk < RestRound; kk++)
        {
            PolyMul(*bigtemp, *temp_mul, S1[90 - kk], S1[91 - kk]);
            PolyAdd(*t1temp, S1[65 - kk], S1[92 - kk]);
            PolyAdd(*t2temp, S2[77 - kk], *t1temp);
            PolyAdd(temps1[kk], *bigtemp, *t2temp);

            PolyMul(*bigtemp, *temp_mul, S2[81 - kk], S2[82 - kk]);
            PolyAdd(*t1temp, S2[68 - kk], S2[83 - kk]);
            PolyAdd(*t2temp, S3[86 - kk], *t1temp);
            PolyAdd(temps2[kk], *bigtemp, *t2temp);

            PolyMul(*bigtemp, *temp_mul, S3[108 - kk], S3[109 - kk]);
            PolyAdd(*t1temp, S3[kk + 45], S3[110 - kk]);
            PolyAdd(*t2temp, S1[68 - kk], *t1temp);
            PolyAdd(temps3[kk], *bigtemp, *t2temp);
        }

        for (int pd = 92; pd >= RestRound; pd--)
            S1[pd].PolyCopy(S1[pd - RestRound]);
        for (int pd = 83; pd >= RestRound; pd--)
            S2[pd].PolyCopy(S2[pd - RestRound]);
        for (int pd = 110; pd >= RestRound; pd--)
            S3[pd].PolyCopy(S3[pd - RestRound]);

        for (int pd = 0; pd < RestRound; pd++)
        {
            S1[pd].PolyCopy(temps3[RestRound - 1 - pd]);
            S2[pd].PolyCopy(temps1[RestRound - 1 - pd]);
            S3[pd].PolyCopy(temps2[RestRound - 1 - pd]);
        }

        cout << "Round : " << SubCount * BlockSize + RestRound << endl;
    }

    delete [] temps1;
    delete [] temps2;
    delete [] temps3;
    delete bigtemp;
    delete temp_mul;
    delete t1temp;
    delete t2temp;
#endif

#if 0 // 写入文件
    clock_t all_start = clock();
    cout << "Write to file now..." << endl;
    int ppt = 0;
    for (int8_t pd = 0; pd < 93; pd++)
    {
        string message = " state_" + to_string(pd);
        S1[pd].write_poly(Round, file, message);
    }
    for (int8_t pd = 0; pd < 84; pd++)
    {
        string message = " state_" + to_string(pd + 93);
        S2[pd].write_poly(Round, file, message);
    }
    for (int8_t pd = 0; pd < 111; pd++)
    {
        string message = " state_" + to_string(pd + 177);
        S3[pd].write_poly(Round, file, message);
    }
    clock_t all_end = clock();
    cout << "Trivium Poly write time : " << (double)(all_end - all_start) / CLOCKS_PER_SEC << endl;
#endif
}
