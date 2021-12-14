// ----------------------------------------------------------------------------------------------------
// 
//              功能: 对给定cube在固定轮数进行分割检测, 保留可能线性或者二次的候选cube, 用于恢复超多项式
//              要求: 要求cube的维数为偶数
// 
// ----------------------------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <ctime>
#include <bitset>
#include <vector>
#include <map>

using namespace std;
typedef uint32_t Type;

// ------------------------------------------------------------------------------------
// Note that the Dim is adjusted according to the given cube and the RetNum = c(n,m).
const Type Dim = 24;  // (Dim = dim(cube)/ 2)
const Type RetNum = 2325;  // (20~23, 2048), (19~22, 1794), (19~21, 232), (20~24, 12951), (19~23, 10903), (21~24, 2325)
Type M[2] = { Dim - 3, Dim };
// ------------------------------------------------------------------------------------

struct IvTerm
{
    bool Term[80] = { 0 };
    void show() const
    {
        for (int i = 0; i < 80; i++)
            if (Term[i])
                cout << i << ' ';
        cout << endl;
    }
};
Type *LocationArr;
typedef bitset<RetNum> TarArr;
typedef map<IvTerm, int16_t> CandidateCubeMap;
bool operator<(const IvTerm &t1, const IvTerm &t2);

// Auxiliary function
Type CnmNumber(Type n, Type m);  // The dimension of a given large cube N, and the dimension of the smallest subcube M to be chosen
Type C_nm(Type arr[], Type n, Type m, Type value = 0, Type deep = 0, Type count = 0);
Type* GetCnmNumner(Type n, Type m[]);
void Moebius_Transform(uint16_t Truth_table[], uint32_t length);

// Partition Test
void GetDiffSplitCube(uint8_t SplitKinds[][2 * Dim], uint32_t SplitTestNum, uint8_t BigCube[], int diffvarnum = 0);
bool Linear_test(CandidateCubeMap &candidatecube, uint16_t RetArr[], uint32_t Ivbits[], uint8_t BigCube[], const int BigDim, uint32_t InitRound, char file[]);
bool Test_Split_SubCube(uint16_t RetArr[], uint32_t Ivbits[], uint8_t BigCube[], const int BigDim, uint32_t InitRound, char file1[], char file2[]);
void SumCube(TarArr ResArr[], uint16_t RetArr[], uint32_t Keybits[], uint32_t IVBits[], uint8_t Cube[], const int dim, uint32_t InitRound);
uint16_t Trivium(uint32_t Keybits[], uint32_t Ivbits[], uint32_t InitRound);
void GetRandomValue(uint32_t Key[], int KEY_num);
void GetRandomCube(uint8_t SubCube1[], uint8_t SubCube2[], const int dim, int count, bool flag = 0);



int main()
{

    uint32_t Round = 815;
    const int TestCubeNum = 1;

    uint8_t BigCube[TestCubeNum][2 * Dim] = {
        { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 15, 16, 17, 18, 19, 22, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36, 38, 40, 42, 43, 46, 51, 53, 55, 57, 60, 61, 62, 64, 66, 68, 70, 72, 78 } 
    };

    uint32_t Ivbits[3] = { 0 };
    uint16_t *RetArr = new uint16_t[uint64_t(1) << Dim];
    LocationArr = GetCnmNumner(Dim, M);

    for (int count = 0; count < TestCubeNum; count++)
    {
        string file_str1 = "PartitionTestResult(Cube_" + to_string(0) + ' ' + to_string(Round) + "r " + to_string(2 * Dim) + "d).txt";
        char file1[50];
        file_str1.copy(file1, file_str1.length(), 0);
        *(file1 + file_str1.length()) = '\0';

        string file_str2 = "DifferentPartitions(Cube_" + to_string(0) + ' ' + to_string(Round) + "r " + to_string(2 * Dim) + "d).txt";
        char file2[50];
        file_str2.copy(file2, file_str2.length(), 0);
        *(file2 + file_str2.length()) = '\0';

        bool PassFlag = Test_Split_SubCube(RetArr, Ivbits, BigCube[count], 2 * Dim, Round, file1, file2);
        cout << "- No." << count << " BigCube have been tested. PassFlag = " << PassFlag << endl;
    }
    
    system("pause");
    return 0;
}


bool Test_Split_SubCube(uint16_t RetArr[], uint32_t Ivbits[], uint8_t BigCube[], const int BigDim, uint32_t InitRound, char file1[], char file2[])
{
    fstream fout;
    CandidateCubeMap candidatecube;

    fout.open(file1, ios_base::app);
    fout << "BigCube(" << BigDim << "): ";
    for (int i = 0; i < BigDim - 1; i++)
        fout << int(BigCube[i]) << ", ";
    fout << int(BigCube[BigDim - 1]) << endl;
    fout.close();

    bool TestRet = Linear_test(candidatecube, RetArr, Ivbits, BigCube, BigDim, InitRound, file2);

    fout.open(file1, ios_base::app);
    fout << "Find Candidate Cubes: " << candidatecube.size() << endl;
    if (TestRet)
    {
        int count = 0;
        for (CandidateCubeMap::iterator pt = candidatecube.begin(); pt != candidatecube.end(); pt++)
        {
            fout << "-No." << count++ << ": ";
            for (int i = 0; i < 80; i++)
                if (pt->first.Term[i])
                    fout << i << ", ";
            fout << "SucceedCount: " << int(pt->second) << endl;
            fout << endl;
        }
    }
    else
        fout << "Find Null" << endl;

    fout << "--------------------------------------------------------------------------------------" << endl << endl;
    fout.close();
    return TestRet;
}



bool Linear_test(CandidateCubeMap &SubCubeSet, uint16_t RetArr[], uint32_t Ivbits[], uint8_t BigCube[], const int BigDim, uint32_t InitRound, char file[])
{
    const int SplitTestNum = 1200;
    const int TestNum = 150;
    const int SmallDim = BigDim / 2;
    uint32_t FailureCount = 0;
    uint32_t KeyArr[3 * (TestNum + 1)] = { 0 };
    TarArr FlagLinearBits1[12], FlagZeroBits1[12], FlagLinearBits2[12], FlagZeroBits2[12];
    TarArr FlagQuadraticBits1[12], FlagOneBits1[12], FlagQuadraticBits2[12], FlagOneBits2[12];
    TarArr output10[12], output20[12], output1[12], output2[12], output3[12], output4[12], output5[12], output6[12], output7[12];
    uint8_t SplitKinds[SplitTestNum][2 * Dim] = { 0 };

    GetDiffSplitCube(SplitKinds, SplitTestNum, BigCube, 4/2);
    fstream fout;
    fout.open(file, ios_base::app);
    for (int i = 0; i < SplitTestNum; i++)
    {
        for (int j = 0; j < 2 * Dim; j++)
            fout << int(SplitKinds[i][j]) << ' ';
        fout << endl;
    }
    fout.close();

    for (int count = 0; count < SplitTestNum; count++)
    {
        bool FailureFlag = 0;
        cout << "Split Cube Count = " << count << endl;
        for (uint8_t i = 0; i < BigDim; i++)
            cout << int(SplitKinds[count][i]) << ' ';
        cout << endl;

        GetRandomValue(KeyArr, TestNum);
        SumCube(output10, RetArr, KeyArr + 3 * TestNum, Ivbits, SplitKinds[count], SmallDim, InitRound);
        SumCube(output20, RetArr, KeyArr + 3 * TestNum, Ivbits, SplitKinds[count] + SmallDim, SmallDim, InitRound);

        for (int i = 0; i < 12; i++)
        {
            FlagLinearBits1[i].reset();
            FlagLinearBits2[i].reset();
            FlagQuadraticBits1[i].reset();
            FlagQuadraticBits2[i].reset();
            FlagZeroBits1[i] = output10[i];
            FlagZeroBits2[i] = output20[i];
            FlagOneBits1[i] = output10[i];
            FlagOneBits2[i] = output20[i];
        }

        for (int pt = 0; pt < TestNum; pt += 3)
        {
            cout << "Linear Test No." << pt / 3;
            uint32_t Key12[3] = { KeyArr[3 * pt] ^ KeyArr[3 * (pt+1)], KeyArr[3 * pt + 1] ^ KeyArr[3 * (pt+1) + 1], KeyArr[3 * pt + 2] ^ KeyArr[3 * (pt+1) + 2] };
            uint32_t Key23[3] = { KeyArr[3 * (pt+1)] ^ KeyArr[3 * (pt+2)], KeyArr[3 * (pt+1) + 1] ^ KeyArr[3 * (pt+2) + 1], KeyArr[3 * (pt+1) + 2] ^ KeyArr[3 * (pt+2) + 2] };
            uint32_t Key13[3] = { KeyArr[3 * pt] ^ KeyArr[3 * (pt + 2)], KeyArr[3 * pt + 1] ^ KeyArr[3 * (pt + 2) + 1], KeyArr[3 * pt + 2] ^ KeyArr[3 * (pt + 2) + 2] };
            uint32_t Key123[3] = { KeyArr[3 * pt] ^ KeyArr[3 * (pt + 1)] ^ KeyArr[3 * (pt + 2)],  KeyArr[3 * pt + 1] ^ KeyArr[3 * (pt + 1) + 1] ^ KeyArr[3 * (pt + 2) + 1],  KeyArr[3 * pt + 2] ^ KeyArr[3 * (pt + 1) + 2] ^ KeyArr[3 * (pt + 2) + 2] };

            SumCube(output1, RetArr, KeyArr + 3 * pt, Ivbits, SplitKinds[count], SmallDim, InitRound);
            SumCube(output2, RetArr, KeyArr + 3 * (pt + 1), Ivbits, SplitKinds[count], SmallDim, InitRound);
            SumCube(output3, RetArr, KeyArr + 3 * (pt + 2), Ivbits, SplitKinds[count], SmallDim, InitRound);
            SumCube(output4, RetArr, Key12, Ivbits, SplitKinds[count], SmallDim, InitRound);
            SumCube(output5, RetArr, Key23, Ivbits, SplitKinds[count], SmallDim, InitRound);
            SumCube(output6, RetArr, Key13, Ivbits, SplitKinds[count], SmallDim, InitRound);
            SumCube(output7, RetArr, Key123, Ivbits, SplitKinds[count], SmallDim, InitRound);

            for (int i = 0; i < 12; i++)
            {
                FlagQuadraticBits1[i] = (output1[i] ^ output2[i] ^ output3[i] ^ output4[i] ^ output5[i] ^ output6[i] ^ output7[i] ^ output10[i]);
                FlagLinearBits1[i] |= (output1[i] ^ output2[i] ^ output4[i] ^ output10[i]);
                FlagLinearBits1[i] |= (output1[i] ^ output3[i] ^ output6[i] ^ output10[i]);
                FlagLinearBits1[i] |= (output2[i] ^ output3[i] ^ output5[i] ^ output10[i]);
                FlagZeroBits1[i] |= output1[i];
                FlagZeroBits1[i] |= output2[i];
                FlagZeroBits1[i] |= output3[i];
                FlagZeroBits1[i] |= output4[i];
                FlagZeroBits1[i] |= output5[i];
                FlagZeroBits1[i] |= output6[i];
                FlagZeroBits1[i] |= output7[i];
                FlagOneBits1[i] &= output1[i];
                FlagOneBits1[i] &= output2[i];
                FlagOneBits1[i] &= output3[i];
                FlagOneBits1[i] &= output4[i];
                FlagOneBits1[i] &= output5[i];
                FlagOneBits1[i] &= output6[i];
                FlagOneBits1[i] &= output7[i];
            }

            SumCube(output1, RetArr, KeyArr + 3 * pt, Ivbits, SplitKinds[count] + SmallDim, SmallDim, InitRound);
            SumCube(output2, RetArr, KeyArr + 3 * (pt + 1), Ivbits, SplitKinds[count] + SmallDim, SmallDim, InitRound);
            SumCube(output3, RetArr, KeyArr + 3 * (pt + 2), Ivbits, SplitKinds[count] + SmallDim, SmallDim, InitRound);
            SumCube(output4, RetArr, Key12, Ivbits, SplitKinds[count] + SmallDim, SmallDim, InitRound);
            SumCube(output5, RetArr, Key23, Ivbits, SplitKinds[count] + SmallDim, SmallDim, InitRound);
            SumCube(output6, RetArr, Key13, Ivbits, SplitKinds[count] + SmallDim, SmallDim, InitRound);
            SumCube(output7, RetArr, Key123, Ivbits, SplitKinds[count] + SmallDim, SmallDim, InitRound);
            for (int i = 0; i < 12; i++)
            {
                FlagQuadraticBits2[i] = (output1[i] ^ output2[i] ^ output3[i] ^ output4[i] ^ output5[i] ^ output6[i] ^ output7[i] ^ output20[i]);
                FlagLinearBits2[i] |= (output1[i] ^ output2[i] ^ output4[i] ^ output20[i]);
                FlagLinearBits2[i] |= (output1[i] ^ output3[i] ^ output6[i] ^ output20[i]);
                FlagLinearBits2[i] |= (output2[i] ^ output3[i] ^ output5[i] ^ output20[i]);
                FlagZeroBits2[i] |= output1[i];
                FlagZeroBits2[i] |= output2[i];
                FlagZeroBits2[i] |= output3[i];
                FlagZeroBits2[i] |= output4[i];
                FlagZeroBits2[i] |= output5[i];
                FlagZeroBits2[i] |= output6[i];
                FlagZeroBits2[i] |= output7[i];
                FlagOneBits2[i] &= output1[i];
                FlagOneBits2[i] &= output2[i];
                FlagOneBits2[i] &= output3[i];
                FlagOneBits2[i] &= output4[i];
                FlagOneBits2[i] &= output5[i];
                FlagOneBits2[i] &= output6[i];
                FlagOneBits2[i] &= output7[i];
            }

            for (int i = 0; i < 6; i++)
                if (((FlagLinearBits1[2 * i].count() == RetNum) && (FlagLinearBits2[2 * i + 1].count() == RetNum)) || ((FlagLinearBits1[2 * i + 1].count() == RetNum) && (FlagLinearBits2[2 * i].count() == RetNum)))
                {
                    FailureCount += 1;
                    cout << "FailureCount = " << FailureCount << endl;
                    FailureFlag = 1;
                    if (FailureCount > 5)
                        return false;
                    break;
                }
            if (FailureFlag)
                break;
            cout << " Pass!" << endl;
        }
        if (FailureFlag)
            continue;
        for (int i = 0; i < 12; i++)
        {
            FlagOneBits1[i].flip();
            FlagOneBits2[i].flip();
        }

        for (int pt = 0; pt < RetNum; pt++)
            for (int pd = 0; pd < RetNum; pd++)
            {
                bool temp = 0;
                for (int k = 0; k < 6; k++)
                {
                    if ((FlagZeroBits1[2 * k].test(pt) & FlagZeroBits2[2 * k + 1].test(pd)) == 0)
                        temp |= 0;
                    else if (((FlagOneBits1[2 * k].test(pt) | FlagQuadraticBits2[2 * k + 1].test(pd)) == 0) || ((FlagQuadraticBits1[2 * k].test(pt) | FlagOneBits2[2 * k + 1].test(pd)) == 0))
                        temp |= 0;
                    else
                        temp |= FlagLinearBits1[2 * k].test(pt) | FlagLinearBits2[2 * k + 1].test(pd);

                    if ((FlagZeroBits1[2 * k + 1].test(pt) & FlagZeroBits2[2 * k].test(pd)) == 0)
                        temp |= 0;
                    else if (((FlagOneBits1[2 * k + 1].test(pt) | FlagQuadraticBits2[2 * k].test(pd)) == 0) || ((FlagQuadraticBits1[2 * k + 1].test(pt) | FlagOneBits2[2 * k].test(pd)) == 0))
                        temp |= 0;
                    else
                        temp |= FlagLinearBits1[2 * k + 1].test(pt) | FlagLinearBits2[2 * k].test(pd);
                }

                // Set the threshold M
                if ((temp == 0) && (count <= 41))
                {
                    IvTerm tempcube;
                    for (int i = 0; i < SmallDim; i++)
                    {
                        if ((LocationArr[pt] >> (SmallDim - 1 - i)) & 1)
                            tempcube.Term[SplitKinds[count][i]] = 1;
                        if ((LocationArr[pd] >> (SmallDim - 1 - i)) & 1)
                            tempcube.Term[SplitKinds[count][SmallDim + i]] = 1;
                    }
                    CandidateCubeMap::iterator iter = SubCubeSet.find(tempcube);
                    if (iter != SubCubeSet.end())
                        iter->second += 1;
                    else
                        SubCubeSet[tempcube] = 1;
                }
                else if ((temp == 0) && (count > 41))
                {
                    IvTerm tempcube;
                    for (int i = 0; i < SmallDim; i++)
                    {
                        if ((LocationArr[pt] >> (SmallDim - 1 - i)) & 1)
                            tempcube.Term[SplitKinds[count][i]] = 1;
                        if ((LocationArr[pd] >> (SmallDim - 1 - i)) & 1)
                            tempcube.Term[SplitKinds[count][SmallDim + i]] = 1;
                    }
                    CandidateCubeMap::iterator iter = SubCubeSet.find(tempcube);
                    if (iter != SubCubeSet.end())
                        iter->second += 1;
                }
            }
        cout << "Find Candidate Cubes :" << SubCubeSet.size() << ", FailureCount = " << FailureCount << ", ";

        if (count >= 45)
            for (CandidateCubeMap::iterator iter = SubCubeSet.begin(); iter != SubCubeSet.end(); )
            {
                if (iter->second < (count - 40 - FailureCount))
                    SubCubeSet.erase(iter++);
                else
                    iter++;
            }
        cout << "After Delete Candidate Cubes :" << SubCubeSet.size() << endl;
    }
    return true;
}


void SumCube(TarArr ResArr[], uint16_t RetArr[], uint32_t Keybits[], uint32_t IVBits[], uint8_t Cube[], const int dim, uint32_t InitRound)
{
    uint32_t CubeFlag[3] = { 0xffffffff, 0xffffffff, 0xffffffff };
    for (int loc = 0; loc < dim; loc++)
        CubeFlag[(Cube[loc] >> 5)] ^= (1 << (31 - (Cube[loc] & 0x1f)));
    uint32_t Ivbits[3] = { IVBits[0], IVBits[1], IVBits[2] };
    uint32_t Len = (uint32_t(1) << (dim));

    for (uint32_t cubevalue = 0; cubevalue < Len; cubevalue++)
    {
        for (int i = 0; i < 3; i++)
            Ivbits[i] &= CubeFlag[i];
        for (int loc = 0; loc < dim; loc++)
            Ivbits[(Cube[loc] >> 5)] ^= (((cubevalue >> (dim - 1 - loc)) & 1) << (31 - (Cube[loc] & 0x1f)));
        RetArr[cubevalue] = Trivium(Keybits, Ivbits, InitRound);
    }

    Moebius_Transform(RetArr, Len);
    for (uint32_t loc = 0; loc < RetNum; loc++)
        for (int bit = 0; bit < 12; bit++)
            ResArr[bit].set(loc, (RetArr[LocationArr[loc]] >> (11 - bit)) & 1);
}


uint16_t Trivium(uint32_t Keybits[], uint32_t Ivbits[], uint32_t InitRound)
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


void GetRandomCube(uint8_t SubCube1[], uint8_t SubCube2[], const int dim, int count, bool flag)
{
    uint8_t *Cube = new uint8_t[2 * dim];
    srand(time(NULL));
    int pt = 0;

    if (!flag)
        Cube[pt++] = (rand() + count) % 80;

        for (; pt < dim; pt++)
            Cube[pt] = SubCube1[pt];

    for (; pt < 2 * dim; pt++)
    {
        Cube[pt] = rand() % 80;
        for (int i = 0; i < pt; i++)
            if (Cube[pt] == Cube[i])
            {
                pt -= 1;
                break;
            }
    }

    if (!flag)
        for (int i = 0; i < dim; i++)
        {
            SubCube1[i] = Cube[i];
            SubCube2[i] = Cube[dim + i];
        }
    else
        for (int i = 0; i < dim; i++)
            SubCube2[i] = Cube[dim + i];

    delete[] Cube;
}


Type* GetCnmNumner(Type n, Type m[])
{
    Type* Len = new Type[m[1] - m[0] + 1];
    Type RetNum = 0;
    for (uint8_t i = m[0]; i <= m[1]; i++)
    {
        Len[i - m[0]] = RetNum;
        RetNum += CnmNumber(n, i);
    }

    static Type* RetArr = new Type[RetNum];
    for (Type pt = 0; pt <= (m[1] - m[0]); pt++)
        C_nm(RetArr + Len[pt], n, m[0] + pt);
    delete[] Len;

    return RetArr;
}


Type C_nm(Type arr[], Type n, Type m, Type value, Type deep, Type count)
{
    Type len = n - m + 1;

    if (deep == 0)
    {
        Type value = 1;
        value <<= (n - 1);
        for (Type i = 0; i < len; i++)
        {
            count = C_nm(arr, n - 1 - i, m - 1, value, deep + 1, count);
            value >>= 1;
        }
    }
    else if (m == 1)
    {
        for (Type i = 0; i < len; i++)
        {
            Type temp = value;
            temp ^= Type(1 << (len - 1 - i));
            arr[count] = temp;
            count += 1;
        }
    }
    else
    {
        for (Type i = 0; i < len; i++)
        {
            Type temp = value;
            temp ^= (Type(1) << (n - 1 - i));
            count = C_nm(arr, n - 1 - i, m - 1, temp, deep + 1, count);
        }
    }
    return count;
}


Type CnmNumber(Type n, Type m)
{
    uint64_t temp = 1;
    m = (n - m > m) ? m : (n - m);
    for (int i = 0; i < m; i++)
        temp *= n - i;
    for (int i = 0; i < m; i++)
        temp /= (i + 1);
    return Type(temp);
}


void Moebius_Transform(uint16_t Truth_table[], uint32_t length)
{
    // Moebius transform for high order bits, using word ops
    for (uint32_t step = 1; step < length; step <<= 1)
        for (uint32_t pos = 0; pos < length; pos += 2 * step)
            for (uint32_t j = 0; j < step; j++)
                Truth_table[step + pos + j] ^= Truth_table[pos + j];
}


bool operator<(const IvTerm &t1, const IvTerm &t2)
{
    for (int i = 0; i < 80; i++)
    {
        if ((!t1.Term[i]) & t2.Term[i])   // t1.test(i) == 0  && t2.test(i) == 1, return true
            return true;
        else if ((t1.Term[i]) & (!t2.Term[i]))  // t1.test(i) == 1  && t2.test(i) == 0, return false
            return false;
    }
    return false;
}


void GetDiffSplitCube(uint8_t SplitKinds[][2 * Dim], uint32_t SplitTestNum, uint8_t BigCube[], int diffvarnum)
{
    const int BigDim = 2 * Dim;
    for (int i = 0; i < BigDim; i++)
        SplitKinds[0][i] = BigCube[i];
    sort(SplitKinds[0], SplitKinds[0] + BigDim);

    for (int count = 1; count < SplitTestNum; count++)
    {
        for (int i = 0; i < BigDim; i++)
            SplitKinds[count][i] = SplitKinds[0][i];
        bool Splitdiff = false;
        while (!Splitdiff)
        {
            random_shuffle(SplitKinds[count], SplitKinds[count] + BigDim);
            sort(SplitKinds[count], SplitKinds[count] + Dim);
            sort(SplitKinds[count] + Dim, SplitKinds[count] + BigDim);

            bool repeatflag = 1;
            for (int cc = 0; cc < count; cc++)
            {
                int diffvarnum1 = 0, diffvarnum2 = 0;
                bool repeatflag1 = 1, repeatflag2 = 1;

                for (int i = 0; i < Dim; i++)
                {
                    if (SplitKinds[count][i] != SplitKinds[cc][i])
                        diffvarnum1 += 1;
                    if (diffvarnum1 > diffvarnum)
                    {
                        repeatflag1 = 0;
                        break;
                    }
                }
                    
                for (int i = 0; i < Dim; i++)
                {
                    if (SplitKinds[count][i] != SplitKinds[cc][Dim + i])
                        diffvarnum2 += 1;
                    if (diffvarnum2 > diffvarnum)
                    {
                        repeatflag2 = 0;
                        break;
                    }
                }
                    
                repeatflag = repeatflag1 | repeatflag2;
                if (repeatflag)
                    break;
            }
            if (!repeatflag)
                Splitdiff = true;
        }
    }
}
