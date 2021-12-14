#include "Poly.h"
#include "Trivium.h"

void RecoverySuperpoly(int16_t OutPutRound, int Cube[], const int dim, int count);

int main()
{
    int16_t OutPutRound = 815;


	int Cube[][48] = {
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 15, 16, 17, 18, 20, 22, 24, 26, 28, 30, 31, 32, 33, 35, 37, 38, 39, 41, 42, 45, 48, 50, 52, 54, 61, 63, 65, 67, 70, 73, 74, 78 }
    };

	int dim[] = { 44 };

	int CubeNum = 1;

    for (int i = 0; i < 1; i++)
        RecoverySuperpoly(OutPutRound, Cube[i], dim[i], i);


    cout << endl << "OK!" << endl;
    system("pause");
    return 0;
}


void RecoverySuperpoly(int16_t OutPutRound, int Cube[], const int dim, int count)
{
    const int16_t TargetRound = 300;
    string file_str = "ExpandingANF" + to_string(TargetRound) + ".txt";
    char ANFFile[30];
    file_str.copy(ANFFile, file_str.length(), 0);
    *(ANFFile + file_str.length() ) = '\0';

    // ����ʽ��д���ļ�
    char FileName[] = "Superpoly815.txt";

#if 0  // �������Ҫ��ʼ�ͻ��ƵĻ������Խ���δ���ע�͵�
    ExpressRecursivelyForTrivium(OutPutRound, TargetRound, ANFFile);
#endif

    uint8_t IvCon[10] = {0};
    Ivterm CubeFlag = 0;
    for (int i = 0; i < dim; i++)
        CubeFlag.set(Cube[i]);
    TermFw CubeTerm;
    for (uint8_t i = 0; i < dim; i++)
        CubeTerm.pterm[Cube[i] >> 5] ^= 1 << (0x1f ^ Cube[i] & 0x1f);
    CubeTerm.deg = dim;

    file_str = "FilteredANF" + to_string(OutPutRound) + '_' + to_string(count) + ".txt";
    char Outfile[30];
    file_str.copy(Outfile, file_str.length(), 0);
    *(Outfile + file_str.length() ) = '\0';

#if 1  // ����Ѿ�ɸ����һ���ˣ����Խ���δ���ע�͵�
    Poly OutPutANF1 = FilterTerm(ANFFile, CubeFlag, OutPutRound, IvCon, Outfile);
#endif

    int16_t InPutRound = 0;
#if 1  // ����Ѿ����й��ڶ���ɸѡ�ˣ����Խ���δ���ע�͵���Ȼ�󲹳�һ���������ʽ�ĺ���
    Poly OutPutANF2 = FilterPartTerm(InPutRound, Outfile, CubeFlag, IvCon);
    OutPutANF2.write_value(InPutRound, Outfile, "Round " + to_string(InPutRound) + " OutPutANF with " + to_string(OutPutANF2.Size) + "terms ");
#endif
    // Poly OutPutANF2(1);
    // ReadANF(InPutRound, OutPutANF2, "FilteredANF815_*.txt");
    
#if 1  // �ָ�������ʽ
    // �жϻ��������£�һ��������м�״̬���������ж��Ƿ��ܹ����Ƶõ���Ӧ�ĳ�����ʽ

    if ( (InPutRound >= 66) || ((InPutRound > 32) && (OutPutANF2.poly[0].deg >= 46)) )
    {
        PolyFw Temp(1);
        string Cubestr = "- No." + to_string(count) + " Cube: ";
        for (int i = 0; i < 80; i++)
            if (CubeFlag.test(i))
                Cubestr += 'v' + to_string(i);

        if (OutPutANF2.Size != 0)
        {
            Cubestr += ": Very Complicated";
            Temp.write_poly(OutPutRound, FileName, Cubestr);
            cout << "Failed!" << endl;
            return;
        }
        else
        {
            Temp.write_poly(OutPutRound, FileName, Cubestr);
            cout << "Failed!" << endl;
            return;
        }
    }
    
    PolyFw* S1 = new PolyFw [93];
    PolyFw* S2 = new PolyFw [84];
    PolyFw* S3 = new PolyFw [111];
    PolyFw* state = new PolyFw [288];
    char TriviumFile[] = "TriviumANF.txt";
    Trivium(S1, S2, S3, CubeFlag, InPutRound, TriviumFile);
    for (int i = 0; i < 93; i++)
        state[i] = S1[i];
    for (int i = 0; i < 84; i++)
        state[93 + i] = S2[i];
    for (int i = 0; i < 111; i++)
        state[177 + i] = S3[i];
    delete [] S1;
    delete [] S2;
    delete [] S3;

#if 0 // �����̬
    for (int i = 0; i < 288; i++)
        state[i].show(InPutRound, "  state" + to_string(i));
#endif

    PolyFw* PolyAllTerm = new PolyFw (200 * SIZEFw);
    PolyFw* PolyMulBig1 = new PolyFw (200 * SIZEFw);
    PolyFw* PolyMulBig2 = new PolyFw (200 * SIZEFw);
    PolyFw* PolyBigTemp = new PolyFw (200 * SIZEFw);

    OutPutANF2.show(InPutRound, "State Poly");
    // ����״̬����������key��iv�ĳ˻���
    for (uint32_t pt = 0; pt < OutPutANF2.Size; pt++)
    {
        TermFw temp;
        PolyMulBig1->poly[0] = temp;
        PolyMulBig1->Size = 1;
        for (uint16_t pd = 0; pd < 288; pd++)
            if ((OutPutANF2.poly[pt].pterm[pd >> 5] >> (0x1f ^ 0x1f & pd)) & 0x1 )
            {
                PolyMul(*PolyMulBig2, *PolyBigTemp, *PolyMulBig1, state[pd]);
                (*PolyMulBig1).PolyCopy(*PolyMulBig2);
            }

        // ֻ����������ʽ����������������ֻ����������������cubeά�����������ȥ��
        (*PolyMulBig2).SuperPolyTerm(CubeFlag.count());
        (*PolyMulBig2).show(InPutRound, "Term_" + to_string(pt));
        PolyAdd(*PolyBigTemp, *PolyAllTerm, *PolyMulBig2);
        (*PolyAllTerm).PolyCopy(*PolyBigTemp);
    }

    // *PolyAllTerm ����������лָ������ĳ�����ʽ�е���
    (*PolyAllTerm).XorTerm(CubeTerm);
    string Cubestr = "- No." + to_string(count) + " Cube: ";
    for (int i = 0; i < 80; i++)
        if (CubeFlag.test(i))
            Cubestr += 'v' + to_string(i);
    (*PolyAllTerm).write_poly(OutPutRound, FileName, Cubestr);

    delete PolyAllTerm;
    delete PolyMulBig1;
    delete PolyMulBig2;
    delete PolyBigTemp;
    delete [] state;
#endif  // �ָ�������ʽ
}