#pragma once
#include <gurobi_c++.h>
#include <bitset>
#include "Poly.h"

#ifndef TRIVIUM_H_
#define TRIVIUM_H_


typedef bitset<80> Ivterm; // 注意使用bitset时, 低位在右边, 高位在左边, 因此bitset[0]访问的是最右边的bit. 
typedef bitset<160> AllVars;
typedef vector<uint32_t> TermCounter;
typedef vector<AllVars> PolyKey;

// --------------------------------------------- 正向计算给定轮数的 ANF ----------------------------------------------
void Trivium(PolyFw S1[], PolyFw S2[], PolyFw S3[], Ivterm &CubeFlag, int Round, char file[]);

// --------------------------------------------- 利用三子集恢复超多项式 ----------------------------------------------
bool CMP(const AllVars &p1, const AllVars &p2);
void TriviumCoreBDPT(GRBModel &Model, GRBVar* Var, int loc[], int round);
bool TriviumEvalBDPT(PolyKey& ResPoly, uint32_t Pterm[], AllVars &Cube, AllVars &Iv_term, int Round);
PolyKey& recover_supperpoly_BDPT(Poly& InPutANF, Ivterm& Cube, int Round, uint8_t iv_con[], char file[]);
void write2file(PolyKey &Superpoly, char file[], AllVars Cube_term);

// ---------------------------------------------- 利用代数次数估计筛项 -----------------------------------------------  = "FilteredANF.txt"
Poly& FilterTerm(char Infile[], Ivterm &CubeFlag, int16_t OutPutRound, uint8_t IvCon[], char Outfile[]);
Poly& FilterPartTerm(int16_t& InPutRound, char file[], Ivterm &CubeFlag, uint8_t IvCon[], uint32_t ConstrTermNum = 5000);

// ----------------------------------------- 针对给定项估计代数次(MILP CBDP) -----------------------------------------
void TriviumCore(GRBModel &Model, GRBVar* Var, uint8_t VarFlag[], int loc[], int round);
bool TriviumEval(Ivterm& Cube, uint32_t Term[], uint8_t fflag[], int Round);

// ------------------------------------------ 针对给定项估计代数次(数值映射) ------------------------------------------
void Numeric_Mapping(int State[], int SState[], Ivterm& Cube, int Round);
int Deg_Mul(int dA[], int dB[], int dC[], int t, int flag);
int min(int arr[], int len);
int max(int arr[], int len);

// -------------------------------------------------- 多项式的回推 --------------------------------------------------
void InvExpressInit(Poly& ZR, Poly& S1, Poly& S94, Poly& S178, Poly& S1S94, Poly& S1S178, Poly& S94S178, Poly& SAll);
void ExpressOneRound(Poly& OutPutANF, Poly& InPutANF, Poly UpdateFunc[], Poly& TempPoly);
void ExpressRecursivelyForTrivium(const int16_t OutPutRound, const int16_t TargetRound, char file[]);
void ReadANF(int16_t &TargetRound, Poly& InPutANF, char file[]);

#endif