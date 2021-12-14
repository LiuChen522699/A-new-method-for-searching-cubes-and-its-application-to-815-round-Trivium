#include <gurobi_c++.h>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <bitset>
#include <ctime>

using namespace std;
typedef bitset<80> IvTerm;

// search for the starting cube
void InitCube(const int dim, uint32_t Ivbits[], uint32_t InitRound);
uint16_t Linear_test(uint32_t Ivbits[], uint8_t Cube[], const int dim, uint32_t InitRound);
uint16_t SumCube(uint32_t Keybits[], uint32_t IVBits[], uint8_t Cube[], const int dim, uint32_t InitRound);
uint16_t Trivium12bits(uint32_t Keybits[], uint32_t Ivbits[], uint32_t InitRound);
void GetRandomValue(uint32_t Key[], int KEY_num);
void GetRandomCube(uint8_t Cube[], const int dim, int count);

// the method to construct potentially good cubes
void TriviumCore(GRBModel &Model, GRBVar* Var, uint8_t VarFlag[], int loc[], int round);
int DegEvalSingleTerm(IvTerm Cube, IvTerm NonCube, int term[], int termlen, int Round);
int DegEval(IvTerm Cube, IvTerm NonCube, int Round, int TargetDeg = 9999);
IvTerm Strategy(IvTerm Cube, IvTerm NonCube, int Round, char file[]);
IvTerm ConstructCube(uint32_t cube[], int dim, uint8_t ivCon[], int Round, char file[]);


int main()
{
#if 1  // search for the starting cube
    uint32_t Ivbits[3] = { 0 };
    uint32_t InitRound = 815;
    const int dim = 16;
    InitCube(dim, Ivbits, InitRound);
#endif 

#if 0 // construct potentially good cubes
    uint8_t Ivcon[10] = {0};
    int Round = 815;
    const int dim = 16;
    const int cubenum = 10;
    uint32_t Cube[cubenum][dim] = {
        {0, 1, 2, 3, 5, 12, 19, 20, 33, 38, 42, 57, 58, 59, 63, 77},
        {0, 1, 2, 4, 12, 21, 23, 29, 35, 39, 43, 46, 48, 51, 62, 71},
        {0, 1, 3, 12, 17, 21, 29, 36, 38, 43, 47, 52, 53, 57, 63, 71},
        {0, 1, 4, 8, 14, 31, 34, 36, 42, 48, 54, 63, 65, 67, 75, 78},
        {0, 1, 4, 9, 11, 15, 16, 22, 29, 35, 41, 44, 45, 62, 74, 76},
        {0, 1, 4, 13, 14, 20, 21, 24, 26, 33, 38, 39, 43, 49, 64, 66},
        {0, 1, 4, 15, 21, 24, 31, 37, 41, 42, 51, 57, 66, 67, 74, 76},
        {0, 1, 5, 6, 9, 10, 12, 19, 21, 22, 23, 25, 26, 42, 52, 53},
        {0, 1, 6, 7, 17, 19, 26, 32, 39, 64, 65, 69, 71, 74, 76, 79},
        {0, 1, 6, 12, 21, 25, 33, 39, 44, 45, 49, 52, 56, 61, 72, 77}
    }

    for (uint32_t pt = 0; pt < cubenum; pt++)
	{
		cout << "- No." << pt << " ExtendCue: ";
		for (int pd = 0; pd < dim; pd++)
			cout << 'v' << Cube[pt][pd] << ", ";
		cout << endl;
		string file_str = "CandidateBigCube(835 " + to_string(pt) + ").txt";
		char OutPutFile[30];
		file_str.copy(OutPutFile, file_str.length(), 0);
		*(OutPutFile + file_str.length()) = '\0';

		IvTerm CubeNew = ConstructCube(Cube[pt], dim, Ivcon, Round, OutPutFile);
	}
#endif 

	system("pause");
	return 0;
}


IvTerm ConstructCube(uint32_t cube[], int dim, uint8_t ivCon[], int Round, char file[])
{
    IvTerm Cube, NonCube;
    for (int i = 0; i < dim; i++)
        Cube.set(cube[i]);
    for (int i = 0; i < 80; i++)
        NonCube.set(i, (ivCon[i>>3] >> (7 ^ i & 7))&1 );
    
    IvTerm NewCube = Strategy(Cube, NonCube, Round, file);
    return NewCube;
}


IvTerm Strategy(IvTerm Cube, IvTerm NonCube, int Round, char file[])
{
    fstream fout;
    int TestNonCube[80][2] = {0};
    IvTerm IvALL = Cube;
    IvALL.flip();
    int TargetDeg = 999;
    int TargetDim = 46;
    bool whilefalg = true;

    while ((Cube.count() < TargetDim) && (TargetDeg > 1) && (whilefalg) )
    {
        int MinDeg = 9999;
        int MinVar = 80;
        
        fout.open(file, ios_base::app);
        fout << "Strategy one!" << endl;
        fout.close();

        for (int pt = 0, count = 0; pt < 80; pt++)
        {
            if (!IvALL.test(pt))
                continue;
            
            TestNonCube[count][1] = pt;
            IvTerm TempCube = Cube;
            TempCube.set(pt, 1);
            cout << "Add Cube Var: " << pt << endl;
            TestNonCube[count][0] = DegEval(TempCube, NonCube, Round, MinDeg);

            if ((TestNonCube[count][0] > 0) && (TestNonCube[count][0] < MinDeg))
            {
                MinDeg = TestNonCube[count][0];
                MinVar = TestNonCube[count][1];
            }
            if (TestNonCube[count][0] <= 0)
                whilefalg = false;
            count += 1;
        }

        fout.open(file, ios_base::app);
        fout << "-------------------- Round: " << Round << ", dim = " << Cube.count() << " --------------------" << endl;
        for (int pt = 0; pt < IvALL.count(); pt++)
        {
            fout << "- Cube: ";
            for (int i = 0; i < 80; i++)
                if (Cube.test(i))
                    fout << i << ' ';
            fout << TestNonCube[pt][1] << ", Deg = " << TestNonCube[pt][0] << endl;
        }
        fout.close();

        Cube.set(MinVar, 1);
        IvALL.reset(MinVar);
        TargetDeg = MinDeg;
    }

    if (TargetDeg == 1)
        return Cube;
    
    fout.open(file, ios_base::app);
    fout << "Strategy two!" << endl;
    fout.close();

    while((Cube.count() < TargetDim) && (TargetDeg > 1) && (whilefalg))
    {
        for (int pt = 0, count = 0; pt < 80; pt++)
        {
            if (!IvALL.test(pt))
                continue;
            
            TestNonCube[count][1] = pt;
            IvTerm TempCube = Cube;
            TempCube.set(pt, 1);
            TestNonCube[count++][0] = DegEval(TempCube, NonCube, Round);
        }

        fout.open(file, ios_base::app);
        fout << "-------------------- Round: " << Round << ", dim = " << Cube.count() << " --------------------" << endl;
        
        int MaxDeg = 0;
        int MaxVar = 80;
        for (int pt = 0; pt < IvALL.count(); pt++)
        {
            fout << "- Cube: ";
            for (int i = 0; i < 80; i++)
                if (Cube.test(i))
                    fout << i << ' ';
            fout << TestNonCube[pt][1] << ", Deg = " << TestNonCube[pt][0] << endl;

            if ((TestNonCube[pt][0] < TargetDeg) && (TestNonCube[pt][0] > MaxDeg) )
            {
                MaxDeg = TestNonCube[pt][0];
                MaxVar = TestNonCube[pt][1];
            }
        }
        fout.close();

        if (MaxVar != 80)
        {
            Cube.set(MaxVar, 1);
            IvALL.reset(MaxVar);
            TargetDeg = MaxDeg;
        }
        else
            break;
    }

    return Cube;
}


int DegEval(IvTerm Cube, IvTerm NonCube, int Round, int TargetDeg)
{
    int Term[120][3] = { {14,38,39}, {13,39,40}, {38,39,40}, {38,39,41}, {59,83,84}, {58,84,85}, {83,84,85}, {83,84,86}, {39,40,118}, {38,39,119}, 
                         {137,149,150}, {136,150,151}, {149,150,151}, {149,150,152}, {84,85,163}, {83,84,164}, {150,151,238}, {149,150,239}, {29,245,246}, {203,245,246}, 
                         {28,246,247}, {202,246,247}, {245,246,247}, {245,246,248}, {44,260,261}, {218,260,261}, {43,261,262}, {217,261,262}, {260,261,262}, {260,261,263}, 
                         {13,14,0}, {25,26,0}, {28,29,0}, {37,38,0}, {14,40,0}, {39,40,0}, {13,41,0}, {43,44,0}, {52,53,0}, {58,59,0}, 
                         {70,71,0}, {59,85,0}, {84,85,0}, {58,86,0}, {106,107,0}, {14,118,0}, {41,118,0}, {13,119,0}, {40,119,0}, {118,119,0}, 
                         {127,128,0}, {133,134,0}, {136,137,0}, {137,151,0}, {150,151,0}, {136,152,0}, {59,163,0}, {86,163,0}, {58,164,0}, {85,164,0}, 
                         {163,164,0}, {172,173,0}, {178,179,0}, {29,202,0}, {28,203,0}, {202,203,0}, {44,217,0}, {43,218,0}, {217,218,0}, {235,236,0}, 
                         {137,238,0}, {152,238,0}, {136,239,0}, {151,239,0}, {238,239,0}, {29,247,0}, {203,247,0}, {246,247,0}, {28,248,0}, {202,248,0}, 
                         {247,248,0}, {44,262,0}, {218,262,0}, {261,262,0}, {43,263,0}, {217,263,0}, {0,0,0}, {12,0,0}, {15,0,0}, {18,0,0}, 
                         {39,0,0}, {42,0,0}, {54,0,0}, {60,0,0}, {72,0,0}, {87,0,0}, {93,0,0}, {105,0,0}, {108,0,0}, {114,0,0}, 
                         {117,0,0}, {129,0,0}, {132,0,0}, {135,0,0}, {138,0,0}, {150,0,0}, {153,0,0}, {159,0,0}, {165,0,0}, {174,0,0}, 
                         {180,0,0}, {192,0,0}, {195,0,0}, {216,0,0}, {219,0,0}, {222,0,0}, {237,0,0}, {240,0,0}, {261,0,0}, {264,0,0}};
    int DimSet[120] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
                       2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
                       2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int degarr[120] = {0};
    int MaxDeg = 0;

    clock_t start = clock();
    for (int i = 0; i < 120; i++)
    {
        degarr[i] = DegEvalSingleTerm(Cube, NonCube, Term[i], DimSet[i], Round);

        cout << "Term" << i << ": ";
        for (int k = 0; k < DimSet[i]; k++)
            cout << 's' << Term[i][k] << ' ';
        cout << ", deg = " << degarr[i] << endl;

        if (degarr[i] > MaxDeg)
            MaxDeg = degarr[i];

        if (degarr[i] >= TargetDeg)
            break;
    }
    clock_t end = clock();
    cout << "Total time of testing a cube: " << double(end - start) / CLOCKS_PER_SEC << 's' << endl;
    
    return MaxDeg;
}



int DegEvalSingleTerm(IvTerm Cube, IvTerm NonCube, int term[], int termlen, int Round)
{
    GRBEnv *env = new GRBEnv();
    GRBModel Model = GRBModel(*env);
    string Model_name = "Cube_Attack_to_Trivium_" + to_string(Round);  
    Model.set(GRB_StringAttr_ModelName, Model_name.c_str()); 
    Model.set(GRB_IntParam_LogToConsole, 0);

    int Loc1[5] = { 65, 170, 90, 91, 92 };
    int Loc2[5] = { 161, 263, 174, 175, 176 };
    int Loc3[5] = { 242, 68, 285, 286, 287 };
    
    const int TargetRound = Round - 200;
    uint8_t iv_flag[80] = { 0 };  
    uint8_t flag[80] = { 0 };  
    
    for (int i = 0; i < 80; i++)
    {
        if (Cube.test(i))
        {
            flag[i] = 1;
            iv_flag[i] = 2;
        }
        else
            iv_flag[i] = NonCube.test(i);
    }

    uint8_t VarFlag[288] = {0}; 
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

    for (int r = 0; r < TargetRound; r++)
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

    GRBLinExpr OutPut; int ZeroFlag = 1;
    for (int i = 0; i < 288; i++)
    {
        bool zeroflag = 1;
        for (int j = 0; j < termlen; j++)
            if (i == term[j])
            {
                ZeroFlag *= VarFlag[i];
                OutPut += s_state[i];
                Model.addConstr(s_state[i] == 1, "final_constr");
                zeroflag = 0;
                break;
            }
        if (zeroflag)
            Model.addConstr(s_state[i] == 0, "fin_round_con_" + to_string(i));
    }
    
    if (!ZeroFlag)
        Model.addConstr(OutPut == 0, "final_constr_Zero");
    
    Model.optimize();
    if (Model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        int ObjVal = int(Model.get(GRB_DoubleAttr_ObjVal));
        delete env;
        return ObjVal;
    }
    else
    {
        delete env;
        return -1;
    }       
}


void TriviumCore(GRBModel &Model, GRBVar* Var, uint8_t VarFlag[], int loc[], int round)
{
    uint8_t temp_Flag[10];
    string var_name[10];
    char var_type[10] = { 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B' };
    double LB[10] = { 0 };
    double UB[10] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

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


void InitCube(const int dim, uint32_t Ivbits[], uint32_t InitRound)
{
	fstream fout;
	
    int Count = 100000;
    uint8_t *cube = new uint8_t[dim];
    uint32_t usefulbits = 0;
	uint32_t pt = 0;

    int cubenum = 1000;
    uint8_t *cubefind = new uint8_t[dim * cubenum];
    int CubeCount = 0;

    string file_str = "FindCube(" + to_string(InitRound) + "r " + to_string(dim) + "d).txt";
    char file[50];
    file_str.copy(file, file_str.length(), 0);
    *(file + file_str.length()) = '\0';

    while (Count && (CubeCount < cubenum))
    {
        GetRandomCube(cube, dim, Count--);
        sort(cube, cube + dim);
        int flag1 = 0;
        for (int cc = 0; cc < CubeCount; cc++)
        {
            int flag2 = 1;
            for (int i = 0; i < dim; i++)
                if (cube[i] != cubefind[CubeCount * dim + i])
                    {flag2 = 0; break;}
            if (flag2)
                {flag1 = 1; break;}
        }
        if (flag1)
            continue;

		cout << "No." << pt++ << ", Cube: ";
		for (int i = 0; i < dim; i++)
			cout << int(cube[i]) << ' ';

        usefulbits = Linear_test(Ivbits, cube, dim, InitRound);
        if (usefulbits)
        {
            for (int i = 0; i < dim; i++)
                cubefind[CubeCount*dim + i] = cube[i];
            CubeCount += 1;
			cout << ": Succeed! Find: " << CubeCount << endl;

			fout.open(file, ios_base::app);
            fout << "{ ";
			for (int i = 0; i < dim-1; i++)
				fout << int(cube[i]) << ", ";
			fout << int(cube[dim-1]) << " }, " << endl;
            fout.close();
        }
		else
			cout << ": Failed!" << endl;
    }

    delete[] cube;
    delete[] cubefind;
}


uint16_t Linear_test(uint32_t Ivbits[], uint8_t Cube[], const int dim, uint32_t InitRound)
{
    uint16_t FlagLinearBits = 0, FlagZeroBits = 0;
    const int TestNum = 100;
    uint32_t *KeyArr = new uint32_t[3 * TestNum];
    GetRandomValue(KeyArr, TestNum);

    uint32_t output0;
    uint32_t KEY_Zreo[3] = { 0 };
    output0 = SumCube(KEY_Zreo, Ivbits, Cube, dim, InitRound);
    FlagZeroBits = output0;

    for (int pt = 0; pt < TestNum; pt += 2)
    {
        uint32_t Key3[3] = { KeyArr[3 * pt] ^ KeyArr[3 * pt + 3], KeyArr[3 * pt + 1] ^ KeyArr[3 * pt + 4], KeyArr[3 * pt + 2] ^ KeyArr[3 * pt + 5] };
        uint16_t output1 = SumCube(KeyArr + 3 * pt, Ivbits, Cube, dim, InitRound);
        uint16_t output2 = SumCube(KeyArr + 3 * (pt + 1), Ivbits, Cube, dim, InitRound);
        uint16_t output3 = SumCube(Key3, Ivbits, Cube, dim, InitRound);
        FlagLinearBits |= (output1 ^ output2 ^ output3 ^ output0);
        FlagZeroBits |= (output1 | output2 | output3);

        if ( FlagLinearBits == 0xFFF)
        {
            delete[] KeyArr;
            return 0;
        }
    }

    delete[] KeyArr;
    return (FlagLinearBits | (FlagZeroBits ^ 0xFFF)) ^ 0xFFF;
}


uint16_t SumCube(uint32_t Keybits[], uint32_t IVBits[], uint8_t Cube[], const int dim, uint32_t InitRound)
{
    uint16_t output = 0;
    uint32_t CubeFlag[3] = { 0xffffffff, 0xffffffff, 0xffffffff };

    for (int loc = 0; loc < dim; loc++)
        CubeFlag[(Cube[loc] >> 5)] ^= (1 << (31 - (Cube[loc] % 32)));
    uint32_t Ivbits[3] = { IVBits[0], IVBits[1], IVBits[2] };

    uint64_t Len = (uint64_t(1) << dim);
    for (uint64_t cubevalue = 0; cubevalue < Len; cubevalue++)
    {
        Ivbits[0] &= CubeFlag[0];
        Ivbits[1] &= CubeFlag[1];
        Ivbits[2] &= CubeFlag[2];
        for (int loc = 0; loc < dim; loc++)
            Ivbits[(Cube[loc] >> 5)] ^= (((cubevalue >> (dim - 1 - loc)) & 1) << (31 - (Cube[loc] % 32)));
        uint16_t temp = Trivium12bits(Keybits, Ivbits, InitRound);
        output ^= temp;
    }
    return output;
}


uint16_t Trivium12bits(uint32_t Keybits[], uint32_t Ivbits[], uint32_t InitRound)
{
    uint32_t State1[3] = { 0 };
    uint32_t State2[3] = { 0 };
    uint32_t State3[4] = { 0, 0, 0, 0xe0000 };
    uint32_t uesbits[12] = { 240, 241, 267, 268, 48, 49, 63, 64, 129, 130, 174, 175 };

    for (int i = 0; i < 3; i++)
    {
        State1[i] = Keybits[i];
        State2[i] = Ivbits[i];
    }
    State1[2] &= 0xFFFF0000;
    State2[2] &= 0xFFFF0000;

    uint32_t RIRs = (InitRound - 111) >> 5;
    uint32_t Shiftbit = (InitRound - 111) & 0x1f;

    for (uint32_t Round = 0; Round < RIRs; Round++)
    {
        uint32_t S171 = (State2[1] << 14) ^ ((State2[2] >> 18) & 0x3fff);
        uint32_t S66 = (State1[1] << 2) ^ ((State1[2] >> 30) & 0x3);
        uint32_t S91 = (State1[1] << 27) ^ ((State1[2] >> 5) & 0x7ffffff);
        uint32_t S92 = (S91 << 1) ^ ((State1[2] & 0x10) > 0);
        uint32_t S93 = (S92 << 1) ^ ((State1[2] & 0x8) > 0);

        uint32_t S162 = (State2[1] << 5) ^ ((State2[2] >> 27) & 0x1f);
        uint32_t S264 = (State3[1] << 23) ^ ((State3[2] >> 9) & 0x7fffff);
        uint32_t S175 = (State2[1] << 18) ^ ((State2[2] >> 14) & 0x3ffff);
        uint32_t S176 = (S175 << 1) ^ ((State2[2] & 0x2000) > 0);
        uint32_t S177 = (S176 << 1) ^ ((State2[2] & 0x1000) > 0);

        uint32_t S69 = (State1[1] << 5) ^ ((State1[2] >> 27) & 0x1f);
        uint32_t S243 = (State3[1] << 2) ^ ((State3[2] >> 30) & 0x3);
        uint32_t S286 = (State3[2] << 13) ^ ((State3[3] >> 19) & 0x1fff);
        uint32_t S287 = (S286 << 1) ^ ((State3[3] & 0x40000) > 0);
        uint32_t S288 = (S287 << 1) ^ ((State3[3] & 0x20000) > 0);

        uint32_t t1 = S66 ^ (S91  & S92) ^ S93  ^ S171;
        uint32_t t2 = S162 ^ (S175 & S176) ^ S177 ^ S264;
        uint32_t t3 = S243 ^ (S286 & S287) ^ S288 ^ S69;

        for (int i = 0; i < 2; i++)
        {
            State1[2 - i] = State1[1 - i];
            State2[2 - i] = State2[1 - i];
            State3[3 - i] = State3[2 - i];
        }
        State3[1] = State3[0];
        State1[0] = t3;
        State2[0] = t1;
        State3[0] = t2;
    }

    uint16_t result = 0;
    for (int i = 0; i < 4; i++)
    {
        result ^= ((State3[(uesbits[i] - Shiftbit - 177) >> 5] >> (0x1f ^ (uesbits[i] - Shiftbit - 177) & 0x1f)) & 0x1) << (11 - i);
        result ^= ((State1[(uesbits[4 + i] - Shiftbit) >> 5] >> (0x1f ^ (uesbits[4 + i] - Shiftbit) & 0x1f)) & 0x1) << (7 - i);
        result ^= ((State2[(uesbits[8 + i] - Shiftbit - 93) >> 5] >> (0x1f ^ (uesbits[8 + i] - Shiftbit - 93) & 0x1f)) & 0x1) << (3 - i);
    }
    return result;
}


void GetRandomValue(uint32_t Key[], int KEY_num)
{
    srand(time(0));

    for (int j = 0; j < KEY_num; j++)
    {
        for (int i = 0; i < 2; i++)
            Key[j * 3 + i] = (rand() << 16) + rand();
        Key[j * 3 + 2] = rand() << 16;
        for (int k = 0; k < j; k++)
        {
            int flag = 1;
            for (int kk = 0; kk < 2; kk++)
                if (Key[j * 3 + kk] != Key[k * 3 + kk])
                {
                    flag = 0; break;
                }
            if (flag)
            {
                j -= 1; break;
            }
        }
    }
}


void GetRandomCube(uint8_t Cube[], const int dim, int count)
{
    srand(time(NULL));
    Cube[0] = (rand() + count) % 80;
    for (int pt = 1; pt < dim; pt++)
    {
        Cube[pt] = rand() % 80;
        for (int i = 0; i < pt; i++)
            if (Cube[pt] == Cube[i])
            {
                pt -= 1;
                break;
            }
    }
}
