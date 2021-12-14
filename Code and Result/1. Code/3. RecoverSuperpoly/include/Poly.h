#pragma once
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>

using namespace std;

const int SIZE = 500000;
const int SIZEFw = 200000;

#ifndef POLY_H_
#define POLY_H_


// ************************************** 回推时状态使用的类，是288个状态组成的项 **************************************
struct Term
{
	uint32_t pterm[9] = {0};
	uint8_t deg = 0;

    void show(char file[], string message)
    {
        fstream fout;
        fout.open(file, ios_base::app);
        fout << "- " << message << ": ";

        if (deg == 0)
            fout << '1' << endl;
        else 
            for (int i = 0; i < 9; i++)
                for (int j = 0; j < 32; j++)
                    if ((pterm[i] >> (31 ^ j)) & 1)
                        fout << 's' << ((i << 5) + j + 1);
        fout << endl;
        fout.close();
    }

    Term& operator=(const Term &p)
    {
        if (this == &p)
            return *this;
        
        for (int i = 0; i < 9; i++)
            this->pterm[i] = p.pterm[i];
        this->deg = p.deg;
        return *this;
    }

};

class Poly
{
public:
	Term* poly;
	uint32_t Size = 0;

	Poly();
	Poly(uint32_t size);
	Poly(const Poly &p);
	~Poly();

	Poly & operator=(const Poly &p);	
	Poly & PolyCopy(const Poly &p);
	Poly & SetPolyLen(const uint32_t len);
    Poly & RemoveDup();

	void write_value(int n, char file[], string message = " ");
	void write_output(int n, char file[], string message = " ");
    void show(int Round, string message = " ");
    
};

Term operator* (const Term &pt1, const Term &pt2);
bool operator< (const Term &p1, const Term &p2);
bool operator<= (const Term &p1, const Term &p2);
bool operator> (const Term &p1, const Term &p2);
bool operator== (const Term &p1, const Term &p2);
bool Term_Greater(const Term &p1, const Term &p2);
bool Divisibility(const Term &BigP, const Term &SmallP);

void PolyAdd(Poly &result, const Poly &p1, const Poly &p2);
void PolyMul(Poly &result, Poly &p1, Term &pt1);
void PolyMul(Poly &result, Poly &Temp_Mul, Poly &p1, Poly &p2);
void PolyMul(Poly &result, Poly &Temp_Mul, Poly& p1, Poly& p2, Poly& p3);


// ************************************** 正向计算时使用的类，是160个Key和Iv组成的项 **************************************
struct TermFw
{
	uint32_t pterm[5] = {0};    // 前80个bit为Iv变元，后80个bit为Key变元
	uint8_t deg = 0;

    void show(char file[], string message)
    {
        fstream fout;
        fout.open(file, ios_base::app);
        fout << "- " << message << ": ";

        if ( (deg == 0) && ( (pterm[0] | pterm[1] | pterm[2] | pterm[3] | pterm[4]) == 0) )
            fout << '1' << endl;
        else 
            for (int i = 0; i < 160; i++)
                if ((pterm[i>>5] >> (0x1f ^ i&0x1f)) & 1)
                {
                    if (i < 80)
                        fout << 'v' << i;
                    else
                        fout << 'k' << i-80; 
                }
        fout << endl;
        fout.close();
    }
};


class PolyFw
{
public:
	TermFw* poly;
	uint32_t Size = 0;

	PolyFw();
	PolyFw(uint32_t size);
	PolyFw(const PolyFw &p);
	~PolyFw();

    PolyFw & SuperPolyTerm(const int deg);
	PolyFw & operator=(const PolyFw &p);	
	PolyFw & PolyCopy(const PolyFw &p);
	PolyFw & SetPolyLen(const uint32_t len);
    PolyFw & XorTerm(const TermFw &term);

    void show(int Round, string message = " "); // 多项式展示到屏幕上
    void write_poly(int n, char file[], string message = " ");
};

TermFw operator* (const TermFw &pt1, const TermFw &pt2);
bool operator< (const TermFw &p1, const TermFw &p2);
bool operator<= (const TermFw &p1, const TermFw &p2);
bool operator> (const TermFw &p1, const TermFw &p2);
bool operator== (const TermFw &p1, const TermFw &p2);
bool Term_Greater(const TermFw &p1, const TermFw &p2);

void PolyAdd(PolyFw &result, const PolyFw &p1, const PolyFw &p2);
void PolyMul(PolyFw &result, PolyFw &p1, TermFw &pt1);
void PolyMul(PolyFw &result, PolyFw &Temp_Mul, PolyFw &p1, PolyFw &p2);
void PolyMul(PolyFw &result, PolyFw &Temp_Mul, PolyFw &p1, PolyFw &p2, PolyFw &p3);






// ************************************** 通用函数，用来估计代数次数 **************************************
uint8_t Weight(uint32_t n);
uint8_t degree(uint32_t pt[], bool flag = 0);
void quickSort(Term s[], int64_t l, int64_t r);
void quickSort(TermFw s[], int64_t l, int64_t r);

#endif