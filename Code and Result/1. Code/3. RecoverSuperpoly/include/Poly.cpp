#include "Poly.h"

// ************************************** ����ʱ״̬ʹ�õ��࣬��288��״̬��ɵ��� **************************************
Poly::Poly()
{
	poly = new Term[SIZE];
	Size = 0;
}

Poly::Poly(uint32_t size)
{
	poly = new Term[size + 1];
	Size = 0;
}

Poly::Poly(const Poly &p)
{
	Size = p.Size;
	poly = new Term[Size + 1];
	for (uint32_t pd = 0; pd < Size; pd++)
		poly[pd] = p.poly[pd];
}

Poly & Poly::operator=(const Poly &p)
{
	if (this == &p)
		return *this;
    delete[] poly;
    Size = p.Size;
    poly = new Term [Size + 1];

    if (Size != 0)
	    for (uint32_t pd = 0; pd < Size; pd++)
		    poly[pd] = p.poly[pd];
	return *this;
}

Poly & Poly::SetPolyLen(const uint32_t len)
{
    delete [] poly;
    poly = new Term [len +1];
    Size = 0;
    return *this;
}

Poly::~Poly()
{
	delete [] poly;
}

// �ø��Ʋ��ı�ԭ������������С������ʹ�õĶ��������ǰ�����
Poly & Poly::PolyCopy(const Poly &p)
{
	if (this == &p)
		return *this;
	Size = p.Size;
	for (uint32_t pd = 0; pd < Size; pd++)
		poly[pd] = p.poly[pd];
	return *this;
}

// ȥ������ʽ���ظ�����
Poly & Poly::RemoveDup()
{
    if (this->Size <= 1)
        return *this;
    
    quickSort(this->poly, int64_t(0), int64_t(this->Size)-1);
    // pd�������е��pf���ǽ�û���ظ�������ǰ��
    uint32_t pd = 0, pf = 0, TempSize = this->Size;
    while (pd < TempSize-1)
    {
        if (this->poly[pd] == this->poly[pd+1])
            pd += 2;
        else
            this->poly[pf++] = this->poly[pd++];
    }
    if (pd == TempSize-1)
        this->poly[pf++] = this->poly[pd];
    this->Size = pf;

    return *this;
}

// ����ʽд���ļ�, ֻ���ÿһ��չ�����������
void Poly::write_output(int Round, char file[], string message)
{
    fstream fout;
    fout.open(file, ios_base::app);

    if (Size == 0)
    {
        fout << "Round_" << Round << " OutPutFunction: " << endl << "# NULL" << message << endl;
        fout.close();
        return ;
    }

    uint32_t length = Size;
    fout << "Round_" << Round << " OutPutFunction with degree = "  
         << int(poly[0].deg) << ", number of terms = "
         << Size << "; " << message << endl << "- Poly = ";
        
    for (uint32_t pr = 0; pr < Size; pr++)
    {
        if (pr)
            fout << " + ";
        for (int i = 0; i < 9; i++)
            for (int j = 0; j < 32; j++)
                if ( (poly[pr].pterm[i] >> (31 ^ j)) & 1 )
                    fout << 's' << ( (i<<5) + j + 1);
    }
    fout << endl << endl;
    fout.close();
}

// ֻ���ÿһ��չ�����������
void Poly::write_value(int Round, char file[], string message)
{
    fstream fout;
    fout.open(file, ios_base::app);

    if (Size == 0)
    {
        fout << "Round_of_Expanding_Recursively_" << Round << " OutPutFunction with degree = " << int(poly[0].deg)
             << ", number of terms = " << Size << "; " << message << endl << "- " << Round << ' ' << Size << endl;
        fout.close();
        return ;
    }

    uint32_t length = Size;
    fout << "Round_of_Expanding_Recursively_" << Round << " OutPutFunction with degree = "  << int(poly[0].deg)
         << ", number of terms = " << Size << "; " << message << endl;
    fout << "- " << Round << ' ' << Size << endl;
        
    for (uint32_t pr = 0; pr < Size; pr++)
    {
        for (int i = 0; i < 9; i++)
            fout << hex << "0x" << poly[pr].pterm[i] << ' ';
        fout << hex << "0x" << int(poly[pr].deg) << endl;
    }
    fout << endl << endl;
    fout.close();
}

// ����ʽ�������Ļ��
void Poly::show(int Round, string message)
{
    if (Size == 0)
    {
        cout << "Round_" << Round << " OutPutFunction: " << endl << "- NULL" << message << endl;
        return ;
    }

    uint32_t length = Size;
    cout << "Round_" << Round << " OutPutFunction with degree = "  
         << int(poly[0].deg) << ", number of terms = "
         << Size << "; " << message << endl << "- Poly = ";
        
    for (uint32_t pr = 0; pr < Size; pr++)
    {
        if (pr)
            cout << " + ";
        for (int i = 0; i < 9; i++)
            for (int j = 0; j < 32; j++)
                if ( (poly[pr].pterm[i] >> (31 ^ j)) & 1 )
                    cout << 's' << ( (i<<5) + j + 1);
    }
    cout << endl;
}

void PolyAdd(Poly &result, const Poly &p1, const Poly &p2)
{
    uint32_t len1 = p1.Size;
    uint32_t len2 = p2.Size;

    if (len1 == 0)
        for (uint32_t pd = 0; pd < len2; pd++)
        {
            result.poly[result.Size] = p2.poly[pd];
            result.Size += 1;
        }
    else if (len2 == 0)
        for (uint32_t pt = 0; pt < len1; pt++)
        {
            result.poly[result.Size] = p1.poly[pt];
            result.Size += 1;
        }
    else if ( (len1 != 0) && (len2 != 0) )
    {
        uint32_t pt = 0; 
        uint32_t pd = 0;
        while( (pt < len1) && (pd < len2) )
        {
            if (p1.poly[pt] > p2.poly[pd])
                result.poly[result.Size++] = p1.poly[pt++];
            else if (p1.poly[pt] < p2.poly[pd])
                result.poly[result.Size++] = p2.poly[pd++];
            else
            {
                pd += 1;
                pt += 1;
            }
        }
        for (; pt < len1; pt++)
            result.poly[result.Size++] = p1.poly[pt];
        for (; pd < len2; pd++)
            result.poly[result.Size++] = p2.poly[pd];
    }
}

void PolyMul(Poly &result, Poly &p1, Term &pt1)
{
    for (uint32_t i = 0; i < p1.Size; i++)
        result.poly[result.Size++] = p1.poly[i] * pt1;
}

// result��Ϊ���ʱ�����ԭ�������ݻᱻ���
void PolyMul(Poly &result, Poly &Temp_Mul, Poly &p1, Poly &p2)
{
    uint32_t len1 = p1.Size;
    uint32_t len2 = p2.Size;
    result.Size = Temp_Mul.Size = 0;

    if ( (len1 != 0) && (len2 != 0) )
        for (uint32_t pd = 0; pd < len2; pd++)
            PolyMul(Temp_Mul, p1, p2.poly[pd]);
    
    // ����ȥ��
    uint32_t temp_len = Temp_Mul.Size;
    if (temp_len > 1)
        quickSort(Temp_Mul.poly, int64_t(0), int64_t(temp_len)-1);

    uint32_t pt = 0;
    if (temp_len > 1)
    {
        while(pt < temp_len-1)
        {
            if (Temp_Mul.poly[pt] == Temp_Mul.poly[pt+1])
                pt += 2;
            else
                result.poly[result.Size++] = Temp_Mul.poly[pt++];
        }
        if (pt == temp_len-1)
            result.poly[result.Size++] = Temp_Mul.poly[pt];
    }
    else if (temp_len == 1)
    {
        result.poly[0] = Temp_Mul.poly[0];
        result.Size = 1;
    }
    Temp_Mul.Size = 0;
}

void PolyMul(Poly &result, Poly &Temp_Mul, Poly& p1, Poly& p2, Poly& p3)
{
    uint32_t TempSize = p1.Size * p2.Size;
    Poly TempPoly(TempSize);
    PolyMul(TempPoly, Temp_Mul, p1, p2);
    PolyMul(result, Temp_Mul, TempPoly, p3);
}

// Term Order and Operation
Term operator* (const Term &pt1, const Term &pt2)
{
	Term result;
    if (pt1.deg == 0)
        result = pt2;
    else if (pt2.deg == 0)
        result = pt1;
    else
    {
        for (int i = 0; i < 9; i++)
		    result.pterm[i] = (pt1.pterm[i]) | (pt2.pterm[i]);
	    result.deg = degree(result.pterm);
    }
	return result;
}

bool operator< (const Term &p1, const Term &p2)
{
	if (p1.deg < p2.deg)
		return true;

	else if (p1.deg > p2.deg)
		return false;

	else
	{
		for (int i = 0; i < 9; i++)
		{
			if (p1.pterm[i] < p2.pterm[i])
				return true;
			else if (p1.pterm[i] > p2.pterm[i])
				return false;
		}
	}
	return false;
}



bool operator<= (const Term &p1, const Term &p2)
{
	if (p1 > p2)
        return false;
    else
        return true;
}

bool operator> (const Term &p1, const Term &p2)
{
	if (p1.deg > p2.deg)
		return true;

	else if (p1.deg < p2.deg)
		return false;

	else
	{
		for (int i = 0; i < 9; i++)
		{
			if (p1.pterm[i] > p2.pterm[i])
				return true;
			else if (p1.pterm[i] < p2.pterm[i])
				return false;
		}
	}
	return false;
}

bool Term_Greater(const Term &p1, const Term &p2)
{
	if (p1.deg > p2.deg)
		return true;

	else if (p1.deg < p2.deg)
		return false;

	else
	{
		for (int i = 0; i < 9; i++)
		{
			if (p1.pterm[i] > p2.pterm[i])
				return true;
			else if (p1.pterm[i] < p2.pterm[i])
				return false;
		}
	}
	return false;
}

bool operator== (const Term &p1, const Term &p2)
{
	// �Ƚϴ���
	if (p1.deg != p2.deg)
		return false;

	// ������ͬʱ���ȽϾ���ֵ
	for (int i = 0; i < 9; i++)
		if (p1.pterm[i] != p2.pterm[i])
			return false;

	// ȫ����ͬ��������
	return true;
}

bool Divisibility(const Term &BigP, const Term &SmallP)
{
    if (BigP.deg < SmallP.deg)
        return false;

    Term temp;
    for (int i = 0; i < 9; i++)
        temp.pterm[i] = BigP.pterm[i] | SmallP.pterm[i];
    if (temp == BigP)
        return true;
    else
        return false;
}


// ************************************** �������ʱʹ�õ��࣬��160��Key��Iv��ɵ��� **************************************
PolyFw::PolyFw()
{
	poly = new TermFw[SIZEFw];
	Size = 0;
}

PolyFw::PolyFw(uint32_t size)
{
	poly = new TermFw[size + 1];
	Size = 0;
}

PolyFw::PolyFw(const PolyFw &p)
{
	Size = p.Size;
	poly = new TermFw[Size + 1];
	for (uint32_t pd = 0; pd < Size; pd++)
		poly[pd] = p.poly[pd];
}

PolyFw & PolyFw::operator=(const PolyFw &p)
{
	if (this == &p)
		return *this;
    delete [] poly;
	Size = p.Size;
	poly = new TermFw [Size + 1];

    if (Size != 0)
        for (uint32_t pd = 0; pd < Size; pd++)
            poly[pd] = p.poly[pd];
	return *this;
}

PolyFw & PolyFw::SetPolyLen(const uint32_t len)
{
    delete [] poly;
    poly = new TermFw [len + 1];
    Size = 0;
    return *this;
}

PolyFw::~PolyFw()
{
	delete [] poly;
}

// �ø��Ʋ��ı�ԭ������������С������ʹ�õĶ��������ǰ�����
PolyFw & PolyFw::PolyCopy(const PolyFw &p)
{
	if (this == &p)
		return *this;
	Size = p.Size;
	for (uint32_t pd = 0; pd < Size; pd++)
		poly[pd] = p.poly[pd];
	return *this;
}

// ����cubeά�������س�����ʽ�е���
PolyFw & PolyFw::SuperPolyTerm(const int deg)
{
    if (this->Size == 0)
        return *this;

    for (uint32_t pd = 0; pd < this->Size; pd++)
        if (this->poly[pd].deg < deg)
        {
            this->Size = pd;
            break;
        }
    return *this;
}

// �����ͬһ���Ŀǰ����ȥ��������ʽ�е�Cube��Ԫ
PolyFw & PolyFw::XorTerm(const TermFw &term)
{
    for (uint32_t pd = 0; pd < this->Size; pd++)
    {
        for (uint8_t pt = 0; pt < 5; pt++)
            this->poly[pd].pterm[pt] ^= term.pterm[pt];
        this->poly[pd].deg = degree(this->poly[pd].pterm, 1);
    }
    return *this;
}

// ����ʽд���ļ�
void PolyFw::write_poly(int Round, char file[], string message)
{
    fstream fout;
    fout.open(file, ios_base::app);

    if (Size == 0)
    {
        fout << "Round_" << Round << " Poly with degree = -1, number of terms = 0; "
             << message << endl << "# NULL" << endl << endl;
        fout.close();
        return ;
    }

    fout << "Round_" << Round << " Poly with degree = " 
            << int(poly[0].deg) << ", number of terms = "
            << Size << "; " << message << endl << "# ";
    for (uint32_t pr = 0; pr < Size; pr++)
    {
        if (pr)
            fout << " + ";
        if ( (pr == Size - 1) && (poly[Size-1].deg == 0) && ( (poly[Size-1].pterm[0]|poly[Size-1].pterm[1]|poly[Size-1].pterm[2]|poly[Size-1].pterm[3]|poly[Size-1].pterm[4]) == 0 ) )
            fout << '1';
        else
            for (int i = 0; i < 160; i++)
                if ( (poly[pr].pterm[i>>5] >> (31 ^ 31&i)) & 1 )
                    {
                        if (i < 80)
                            fout << 'v' << i;
                        else 
                            fout << 'k' << i - 80;    
                    }    
    }
    fout << endl << endl;
    fout.close();
}

// ����ʽչʾ����Ļ��
void PolyFw::show(int Round, string message)
{
    if (Size == 0)
    {
        cout << "Round_" << Round << " Poly with degree = -1, number of terms = 0 : "
             << message << endl << "# NULL" << endl;
        return ;
    }

    cout << "Round_" << Round << " Poly with degree = " 
            << int(poly[0].deg) << ", number of terms = "
            << Size << "; " << message << endl << "- Poly = ";
    for (uint32_t pr = 0; pr < Size; pr++)
    {
        if (pr)
            cout << " + ";
        if ( (pr == Size - 1) && (poly[Size-1].deg == 0) && ( (poly[Size-1].pterm[0]|poly[Size-1].pterm[1]|poly[Size-1].pterm[2]|poly[Size-1].pterm[3]|poly[Size-1].pterm[4]) == 0 ) )
            cout << '1';
        else
            for (int i = 0; i < 160; i++)
                if ( (poly[pr].pterm[i>>5] >> (31 ^ 31&i)) & 1 )
                    {
                        if (i < 80)
                            cout << 'v' << i;
                        else 
                            cout << 'k' << i - 80;    
                    }    
    }
    cout << endl;
}

void PolyAdd(PolyFw &result, const PolyFw &p1, const PolyFw &p2)
{
    uint32_t len1 = p1.Size, len2 = p2.Size;
    result.Size = 0;

    if (len1 == 0)
        for (uint32_t pd = 0; pd < len2; pd++)
            result.poly[result.Size++] = p2.poly[pd];
    else if (len2 == 0)
        for (uint32_t pt = 0; pt < len1; pt++)
            result.poly[result.Size++] = p1.poly[pt];
    else if ( (len1 != 0) && (len2 != 0) )
    {
        uint32_t pt = 0, pd = 0; 
        while( (pt < len1) && (pd < len2) )
        {
            if (p1.poly[pt] > p2.poly[pd])
                result.poly[result.Size++] = p1.poly[pt++];
            else if (p1.poly[pt] < p2.poly[pd])
                result.poly[result.Size++] = p2.poly[pd++];
            else
            {
                pd += 1;
                pt += 1;
            }
        }
        for (; pt < len1; pt++)
            result.poly[result.Size++] = p1.poly[pt];
        for (; pd < len2; pd++)
            result.poly[result.Size++] = p2.poly[pd];
    }
}

void PolyMul(PolyFw &result, PolyFw &p1, TermFw &pt1)
{
    for (uint32_t i = 0; i < p1.Size; i++)
        result.poly[result.Size++] = p1.poly[i] * pt1;
}

void PolyMul(PolyFw &result, PolyFw &Temp_Mul, PolyFw &p1, PolyFw &p2)
{
    uint32_t len1 = p1.Size;
    uint32_t len2 = p2.Size;
    Temp_Mul.Size = result.Size = 0;

    if ( (len1 != 0) && (len2 != 0) )
        for (uint32_t pd = 0; pd < len2; pd++)
            PolyMul(Temp_Mul, p1, p2.poly[pd]);

    // ����ȥ��
    uint32_t temp_len = Temp_Mul.Size;
    if (temp_len > 1)
        quickSort(Temp_Mul.poly, int64_t(0), int64_t(temp_len)-1);

    uint32_t pt = 0;
    if (temp_len > 1)
    {
        while(pt < temp_len-1)
        {
            if (Temp_Mul.poly[pt] == Temp_Mul.poly[pt+1])
                pt += 2;
            else
                result.poly[result.Size++] = Temp_Mul.poly[pt++];
        }
        if (pt == temp_len-1)
            result.poly[result.Size++] = Temp_Mul.poly[pt];
    }
    else if (temp_len == 1)
    {
        result.poly[0] = Temp_Mul.poly[0];
        result.Size = 1;
    }
    Temp_Mul.Size = 0;
}

void PolyMul(PolyFw &result, PolyFw &Temp_Mul, PolyFw& p1, PolyFw& p2, PolyFw& p3)
{
    uint32_t TempSize = p1.Size * p2.Size;
    PolyFw TempPoly(TempSize);
    PolyMul(TempPoly, Temp_Mul, p1, p2);
    PolyMul(result, Temp_Mul, TempPoly, p3);
}

// Term Order and Operation
TermFw operator* (const TermFw &pt1, const TermFw &pt2)
{
	TermFw result;
	for (int i = 0; i < 5; i++)
		result.pterm[i] = (pt1.pterm[i]) | (pt2.pterm[i]);
	result.deg = degree(result.pterm, true);
	return result;
}

bool operator< (const TermFw &p1, const TermFw &p2)
{
	if (p1.deg < p2.deg)
		return true;

	else if (p1.deg > p2.deg)
		return false;

	else
	{
		for (int i = 0; i < 5; i++)
		{
			if (p1.pterm[i] < p2.pterm[i])
				return true;
			else if (p1.pterm[i] > p2.pterm[i])
				return false;
		}
	}
	return false;
}

bool operator<= (const TermFw &p1, const TermFw &p2)
{
	if (p1 > p2)
        return false;
    else
        return true;
}

bool operator> (const TermFw &p1, const TermFw &p2)
{
	if (p1.deg > p2.deg)
		return true;

	else if (p1.deg < p2.deg)
		return false;

	else
	{
		for (int i = 0; i < 5; i++)
		{
			if (p1.pterm[i] > p2.pterm[i])
				return true;
			else if (p1.pterm[i] < p2.pterm[i])
				return false;
		}
	}
	return false;
}

bool Term_Greater(const TermFw &p1, const TermFw &p2)
{
	if (p1.deg > p2.deg)
		return true;

	else if (p1.deg < p2.deg)
		return false;

	else
	{
		for (int i = 0; i < 5; i++)
		{
			if (p1.pterm[i] > p2.pterm[i])
				return true;
			else if (p1.pterm[i] < p2.pterm[i])
				return false;
		}
	}
	return false;
}

bool operator== (const TermFw &p1, const TermFw &p2)
{
	// �Ƚϴ���
	if (p1.deg != p2.deg)
		return false;

	// ������ͬʱ���ȽϾ���ֵ
	for (int i = 0; i < 5; i++)
		if (p1.pterm[i] != p2.pterm[i])
			return false;

	// ȫ����ͬ��������
	return true;
}

// *************************************** ͨ�ú��� ***************************************

void quickSort(Term s[], int64_t l, int64_t r)
{
	if (l < r)
	{      
		int64_t i = l, j = r;
        Term x = s[l];
		while (i < j)
		{
			while((i < j) && (s[j] <= x)) // ���������ҵ�һ������x����
				j--; 
			if(i < j)
				s[i++] = s[j];
			while(i < j && s[i] > x) // ���������ҵ�һ��С�ڵ���x����
				i++; 
			if(i < j)
				s[j--] = s[i];
		}
		s[i] = x;
		quickSort(s, l, i - 1); // �ݹ����
		quickSort(s, i + 1, r);
	}
}

void quickSort(TermFw s[], int64_t l, int64_t r)
{
	if (l < r)
	{      
		int64_t i = l, j = r;
        TermFw x = s[l];
		while (i < j)
		{
			while((i < j) && (s[j] <= x)) // ���������ҵ�һ������x����
				j--; 
			if(i < j)
				s[i++] = s[j];
			while(i < j && s[i] > x) // ���������ҵ�һ��С�ڵ���x����
				i++; 
			if(i < j)
				s[j--] = s[i];
		}
		s[i] = x;
		quickSort(s, l, i - 1); // �ݹ����
		quickSort(s, i + 1, r);
	}
}

// ��ȡ��������,������ȡ��Ĵ���
uint8_t Weight(uint32_t n)
{
	n = n - ((n >> 1) & 0x55555555);
	n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
	n = (n & 0x0f0f0f0f) + ((n >> 4) & 0x0f0f0f0f);
	n = n + (n >> 8);
	n = n + (n >> 16);
	return uint8_t(n & 0x0000003f);
}

// ����Ĵ���
uint8_t degree(uint32_t pt[], bool flag)
{
    uint8_t value = 0;
    if (!flag)  // flag = 0: ��ʾ����ʱʹ�õĴ�������
        for (int i = 0; i < 9; i++)
            value += Weight(pt[i]);
    else    // flag = 1: ��ʾ�������ʱʹ�õĴ�������, ��ʱֻ��cube��Ԫ��Ϊ�����������
        value = Weight(pt[0]) + Weight(pt[1]) + Weight(pt[2] & 0xffff0000);
    
    return value;
}