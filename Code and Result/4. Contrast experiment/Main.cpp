#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <gurobi_c++.h>
#include <bitset>
#include <iomanip>
#include <ctime>

using namespace std;

typedef bitset<160> Term;
typedef vector<uint32_t> TermCounter;
typedef vector<Term> Poly;

bool CMP(const Term &p1, const Term &p2);
void GetRandomCube(uint8_t Cube[], int dim);
void TriviumCore(GRBModel &Model, GRBVar* Var, int loc[], int round);
void TriviumEval(GRBModel &Model, Term &Cube_term, Term &Iv_term, int Round);
int ToSolveTheModel(GRBModel &Model, Poly &Superpoly, TermCounter &termcounter, Term &cube);
void write2file(GRBModel &Model, Poly &Superpoly, TermCounter &termcounter, char file[], double usetime);
void recover_supperpoly(uint8_t Cube[], const int dim, int StartRound, int EndRound, uint8_t iv_con[], int count);


int main()
{
    int StartRound = 815;
    int EndRound = 816;
    uint8_t iv_con[10] = { 0 };

    for (int count = 0; count < 1000; count++)
    {
        uint8_t choseflag = (rand() % 5);
        uint8_t dim = ((choseflag == 0) * 42) + ((choseflag == 1) * 43) + ((choseflag == 2) * 44) + ((choseflag == 3) * 45) + ((choseflag == 4) * 46);
        cout << "dim = " << int(dim) << ", ";

        uint8_t *Cube = new uint8_t[dim];
        GetRandomCube(Cube, dim);

        for (int j = 0; j < dim; j++)
            cout << int(Cube[j]) << ' ';
        cout << endl;

        recover_supperpoly(Cube, dim, StartRound, EndRound, iv_con, count);
    }

    system("pause");
    return 0;
}


void recover_supperpoly(uint8_t Cube[], const int dim, int StartRound, int EndRound, uint8_t iv_con[], int count)
{
    char file[] = "result(815).txt";
    string log_file = "log_file_gurobi.log";
    GRBEnv env = GRBEnv(log_file);

	// do a systematic search for the k- best solutions
    env.set(GRB_IntParam_PoolSearchMode, 2);
	// Limit how many solutions to collect
    env.set(GRB_IntParam_PoolSolutions, 2000000000);
	// Limit the search space by setting a gap for the worst possible solution that will be acce
    env.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);
    env.set(GRB_IntParam_Threads, 11);
    env.set(GRB_DoubleParam_TimeLimit, 3600*12);

    Term Cube_term(0), Iv_term(0);

	for (int i = 0; i < 10; i++)
        for (int j = 0; j < 8; j++)
            if ( !((iv_con[i] >> (7-j))&0x1) )
                Iv_term.set(80 + 8*i+j, 1);
	
    for (int i = 0; i < dim; i++)
    {
        Cube_term.set(80 + Cube[i], 1);
        Iv_term.reset(80 + Cube[i]);
    }

    fstream fout;
    fout.open(file, ios_base::app);
    fout << "-No." << count << ", Cube(" << int(dim) << "): ";
    for (int i = 0; i < 80; i++)
        if (Cube_term.test(80 + i))
            fout << 'v' << i << ", ";
    fout << endl;
    fout.close();

    for (int round = StartRound; round < EndRound; round++)
    {
        GRBModel Model = GRBModel(env);
        string Model_name = "Cube_Attack_to_trivium_" + to_string(round) + " Round";
		Model.set(GRB_StringAttr_ModelName, Model_name);

        TriviumEval(Model, Cube_term, Iv_term, round);
        
        Poly Superpoly;
        TermCounter termcounter;

        clock_t starttime = clock();
        int flag = ToSolveTheModel(Model, Superpoly, termcounter, Cube_term);
        clock_t endtime = clock();
        cout << "Model Optimize time: " << double((endtime - starttime) / CLOCKS_PER_SEC) << 's' << endl;;
        
        if (flag)
            write2file(Model, Superpoly, termcounter, file, double(endtime - starttime));
        else
        {
            fout.open(file, ios_base::app);
            fout << Model.get(GRB_StringAttr_ModelName)<< endl << "Model Optimize time: " << double(endtime - starttime)
                << "\nTimeLimit\n\n********************************************************\n\n";
            fout.close();
        }
    }
}


void write2file(GRBModel &Model, Poly &Superpoly, TermCounter &termcounter, char file[], double usetime)
{
    fstream fout;
    fout.open(file, ios_base::app);
    fout << Model.get(GRB_StringAttr_ModelName) << endl
         << "Model Optimize time: " << usetime << "\nSuperpoly = ";

    int term_flag = 0;
    if (Superpoly.size())
    {
        for (int loc = 0; loc < Superpoly.size(); loc++)
            if (termcounter[loc]&0x1)
            {
                if (term_flag)
                    fout << " + ";
                term_flag += 1;

                if (Superpoly[loc] == 0)
                    fout << '1';
                else
                    for (int i = 0; i < 160; i++)
                        if (Superpoly[loc][i] && (i < 80))
                            fout << 'x' << i+1;
                        else if (Superpoly[loc][i] && (i >= 80))
                            fout << 'v' << i-79;
            }
        
        if (!term_flag)
            fout << 0;
        
        fout << endl << "Superpoly Size: " << term_flag << "\n------------------------------------------------\n";
        for (int loc = 0; loc < Superpoly.size(); loc++)
        {
            fout << "- " << setw(8) << termcounter[loc];
            fout << "  |    " << (termcounter[loc] & 0x1) << "    |  ";
			if (Superpoly[loc] == 0)
				fout << 1 << endl;
			else
			{
				for (int i = 0; i < 160; i++)
					if (Superpoly[loc][i] && (i < 80))
						fout << 'x' << i+1;
					else if (Superpoly[loc][i] && (i >= 80))
						fout << 'v' << i-79;
				fout << endl;
			}
        }

        fout << "------------------------------------------------\n\n********************************************************\n\n";
        fout.close();
    }
    else
        fout << 0 << endl << "Superpoly Size: " << term_flag << "\n\n********************************************************\n\n";
}


// 返回0表示在限制时间内没有求解完成；返回1表示有解；返回2表示无解
int ToSolveTheModel(GRBModel &Model, Poly &Superpoly, TermCounter &termcounter, Term &cube)
{
    fstream fout;
    Poly temppoly;
    Superpoly.clear();
    termcounter.clear();

    clock_t starttime = clock();
    Model.optimize();
    clock_t endtime = clock();
    cout << "Model Optimize time: " << double((endtime - starttime)/CLOCKS_PER_SEC) << 's' << endl;;

    if (Model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
    {
        cout << "TimeLimit" << endl;
        return 0;
    }

    // 在有解的情况下, 输出所有的解
    if (Model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        int nSolutions = Model.get(GRB_IntAttr_SolCount);
        cout << " Number of solutions found : " << nSolutions << endl;

        for (int res = 0; res < nSolutions; res++)
        {
            Term p_term(0);
            Model.set(GRB_IntParam_SolutionNumber, res);

            for (int loc = 0; loc < 80; loc++)
            {
                GRBVar TempVar1 = Model.getVarByName("s_init_" + to_string(loc));
                GRBVar TempVar2 = Model.getVarByName("s_init_" + to_string(loc + 93));

                if (TempVar1.get(GRB_DoubleAttr_Xn) == 1)
                    p_term.set(loc);
                if (TempVar2.get(GRB_DoubleAttr_Xn) == 1)
                    p_term.set(80 + loc);
            }
            p_term ^= cube;
            temppoly.push_back(p_term);
        }

        // 处理超多项式的重复问题
        sort(temppoly.begin(), temppoly.end(), CMP);

        Superpoly.push_back(temppoly[0]);
        termcounter.push_back(1);
        for (int loc = 1; loc < temppoly.size(); loc++)
        {
            if (*(Superpoly.end() - 1) == temppoly[loc])
                *(termcounter.end() - 1) += 1;
            else
            {
                Superpoly.push_back(temppoly[loc]);
                termcounter.push_back(1);
            }
        }
        temppoly.clear();
        return 1;
    }
    else
        return 2;
}


void TriviumEval(GRBModel &Model, Term &Cube_term, Term &Iv_term, int Round)
{
    int Loc1[5] = { 65, 170, 90, 91, 92 };
    int Loc2[5] = { 161, 263, 174, 175, 176 };
    int Loc3[5] = { 242, 68, 285, 286, 287 };

    GRBVar* s_state;
    string s_init_str[288];
    char s_type[288];
    double LB[288];
    double UB[288];
    for (int i = 0; i < 288; i++)
    {
        s_init_str[i] = "s_init_" + to_string(i);
        s_type[i] = 'B';
        LB[i] = 0;
        UB[i] = 1;
    }
    s_state = Model.addVars(LB, UB, LB, s_type, s_init_str, 288);
   
    GRBLinExpr key_sum;
	for (int i = 0; i < 80; i++)
		key_sum += s_state[i];
    Model.setObjective(key_sum, GRB_MAXIMIZE);

    for (int i = 0; i < 80; i++)
    {
        string Iv_constr = "Iv_constr_" + to_string(i);
        if (Cube_term.test(80+i))
            Model.addConstr(s_state[93 + i] == 1, Iv_constr);
		else if (Iv_term.test(80+i))
            Model.addConstr(s_state[93 + i] == 0, Iv_constr);
    }

    for (int i = 80; i < 93; i++)
        Model.addConstr(s_state[i] == 0, "constant_constr_" + to_string(i));
    for (int i = 173; i < 285; i++)
        Model.addConstr(s_state[i] == 0, "constant_constr_" + to_string(i));

    for (int r = 0; r < Round; r++)
    {
        TriviumCore(Model, s_state, Loc1, r);
        TriviumCore(Model, s_state, Loc2, r);
        TriviumCore(Model, s_state, Loc3, r);
        GRBVar temp = s_state[287];
        for (int i = 287; i > 0; i--)
            s_state[i] = s_state[i-1];
        s_state[0] = temp;
    }

    for (int i = 0; i < 288; i++)
    {
        if ((i == 65) | (i == 92) | (i == 161) | (i == 176) | (i == 242) | (i == 287))
            continue;
        else
            Model.addConstr(s_state[i] == 0, "fin_round_con_" + to_string(i));
    }
    
    Model.addConstr(s_state[65] + s_state[92] + s_state[161] + s_state[176] + s_state[242] + s_state[287] == 1, "final_constr");
}


void TriviumCore(GRBModel &Model, GRBVar* Var, int loc[], int round)
{
    GRBVar* var_new;
    string var_name[10];
    char var_type[10] = {'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B'};
    double UB[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
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
        Model.addConstr(Var[loc[i]] <= var_new[i] + var_new[5 + i] , "Copy1 Constr_" + to_string(round) + '_' + to_string(i));
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


bool CMP(const Term &p1, const Term &p2)
{
    for (int i = 0; i < 160; i++)
        if (p1[i] > p2[i])
            return true;
        else if (p1[i] < p2[i])
            return false;
	return true;
}


void GetRandomCube(uint8_t Cube[], int dim)
{
    srand(time(NULL));

    for (int pt = 0; pt < dim; pt++)
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