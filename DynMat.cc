#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include "DynMat.h"

#if MAT_DEBUG == 1
double& Mat::operator()(int i, int j)
{
    if (i < 0 || i >= nrow) {
        printf(" out of range1 %d %d\n", i, nrow);
        exit(995);
    }
    if (j < 0 || j >= ncol) {
        printf(" out of range2 %d %d\n", j, ncol);
        exit(996);
    }
    return val[i][j];
}

double& Mat::operator ()(int i)
{
    if (i < 0 || i >= nrow) {
        printf(" out of range3 %d %d\n", i, nrow);
        exit(997);
    }
    return val[i][0];
}

double*& Mat::operator [](int i)
{
    if (i < 0 || i >= nrow) {
        printf(" out of range4 %d %d\n", i, nrow);
        exit(998);
    }
    return val[i];
}
#endif

double sq(double x)
{
    return x * x;
};

Mat::Mat(M_TYPE mt, int nrow_, int ncol_, char cLable_[IdentLength])
{
    myType = mt;
    nrow = maxrow = nrow_;
    ncol = maxcol = ncol_;
    strcpy(cLable, cLable_);
    dEpsilon = 1e-8;
    dEpsilonReg = 2e-8;
    if (mt == HEAP) {
        storage = new double[nrow * ncol];
        row = new bool[nrow];
        col = new bool[ncol];
    }
    if (mt == MALLOC) {
        storage = (double*)malloc(sizeof(double) * nrow * ncol);
        row = (bool*)malloc(sizeof(bool) * nrow);
        col = (bool*)malloc(sizeof(bool) * ncol);
    }
    for (register int i = 0; i < nrow; i++) row[i] = true;
    for (register int i = 0; i < ncol; i++) col[i] = true;
    if (mt == HEAP)
        val = new double* [nrow];
    if (mt == MALLOC)
        val = (double**)malloc(sizeof(double*) * nrow);
    for (register int i = 0; i < nrow; i++) val[i] = storage + i * ncol;
    myNum_ = -1;
}

Mat::Mat(int nrow_, int ncol_, char cLable_[IdentLength])
{
    myType = HEAP;
    nrow = maxrow = nrow_;
    ncol = maxcol = ncol_;
    strcpy(cLable, cLable_);
    dEpsilon = 1e-8;
    dEpsilonReg = 2e-8;
    storage = new double[nrow * ncol];
    row = new bool[nrow];
    col = new bool[ncol];

    for (register int i = 0; i < nrow; i++) row[i] = true;
    for (register int i = 0; i < ncol; i++) col[i] = true;
    val = new double* [nrow];

    for (register int i = 0; i < nrow; i++) val[i] = storage + i * ncol;
    myNum_ = -1;
}

Mat::Mat(const Mat& M)
{
    nrow = M.nrow; ncol = M.ncol;
    maxrow = M.maxrow; maxcol = M.maxcol;
    dEpsilon = M.dEpsilon;
    dEpsilonReg = M.dEpsilonReg;
    strcpy(cLable, M.cLable);
    storage = new double[nrow * ncol];
    memcpy((void*)storage, (void*)M.storage, sizeof(double) * nrow * ncol);
    row = new bool[nrow];
    col = new bool[ncol];
    for (register int i = 0; i < nrow; i++) row[i] = M.row[i];
    for (register int i = 0; i < ncol; i++) col[i] = M.col[i];
    val = new double* [nrow];
    for (register int i = 0; i < nrow; i++) val[i] = storage + i * ncol;
};

Mat::~Mat()
{
    //  printf("%d\n", myType);
    if (myType == HEAP) {
        delete[]val;
        delete[]storage;
        delete[]row;
        delete[]col;
    }
    if (myType == MALLOC) {
        if (val)
            free(val);
        if (storage)
            free(storage);
        if (row)
            free(row);
        if (col)
            free(col);
    }
}

int Mat::nzRow()
{
    int nz = 0;
    for (register int i = 0; i < nrow; i++)
        if (row[i]) nz++;
    return nz;
}

void Mat::Reset(int nrow_, int ncol_)
{
    nrow = nrow_; ncol = ncol_;
    //  if(nrow_ > maxrow){
    //    printf("row %d %d\n", nrow_, maxrow);
    //    exit(512);
    //  }
    //  if(ncol_ > maxcol){
    //    printf("col %d %d\n", ncol_, maxcol);
    //    exit(513);
    //  }
    for (register int i = 0; i < nrow; i++) row[i] = true;
    for (register int i = 0; i < ncol; i++) col[i] = true;
    for (register int i = 0; i < nrow; i++) val[i] = storage + i * ncol;
}

void Mat::ResetAvail(bool br, bool bc)
{
    for (register int i = 0; i < nrow; i++) row[i] = br;
    for (register int i = 0; i < ncol; i++) col[i] = bc;
}

void Mat::Move(const Mat& M)
//void Mat::Move(Mat & M)
{
    Mat& S = *this;
    for (register int i = 0; i < M.nrow; i++)
        S.row[i] = M.row[i];
    //  memcpy( (void*)S.storage, (void*)M.storage, sizeof(double) * M.nrow * M.ncol );
    for (register int i = 0; i < M.ncol; i++)
        S.col[i] = M.col[i];
    for (register int i = 0; i < M.nrow; i++)
        for (register int j = 0; j < M.ncol; j++)
            S(i, j) = M.val[i][j];// M(i, j);

      //  strcpy(cLable, M.cLable);
};

void Mat::SetSub(int ib, int jb, Mat& M,
    int imb, int ime, int jmb, int jme)
{
    Mat& S = *this;
    for (register int i = imb; i <= ime; i++)
        for (register int j = jmb; j <= jme; j++)
            S(ib + i - imb, jb + j - jmb) = M(i, j);
};

void Mat::Zero_()
{
    Mat& S = (*this);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++)
            S(i, j) = 0.;
    }
}
void Mat::Zero()
{
    Mat& S = (*this);
    S.Zero(0, nrow - 1, 0, ncol - 1);
    //  for(int i=0; i < nrow; i++)
    //  {
    //    if(row[i])
    //    for(int j=0; j < ncol; j++)
    //      if(col[j])
    //        S(i,j)=0.;
    //  }
}
void Mat::Zero(int ib, int ie, int jb, int je)
{
    Mat& S = (*this);
    for (register int i = ib; i <= ie; i++)
    {
        if (row[i])
            for (register int j = jb; j <= je; j++)
                if (col[j])
                    S(i, j) = 0.;
    }
};

void Mat::Ident()
{
    Mat& S = (*this);
    S.Zero();
    int n = nrow;
    if (ncol < n)
        n = ncol;
    for (register int i = 0; i < n; i++)
        if (row[i])
            S(i, i) = 1.;
}

double Mat::Scal(Mat& A, int n)
{
    Mat& S = (*this);
    double sc = 0.;
    for (register int i = 0; i < n; i++)
    {
        if (row[i] && A.row[i])
            sc += S(i, 0) * A(i, 0);
    }
    return sc;
}

double Mat::ScalRow(Mat& A, int n)
{
    Mat& S = (*this);
    double sc = 0.;
    for (register int i = 0; i < n; i++)
    {
        if (col[i] && A.col[i])
            sc += S(0, i) * A(0, i);
    }
    return sc;
}

void Mat::ScaleMat(double* w)
{
    ScaleMat(w, nrow);
}
void Mat::ScaleMat2(double* w)
{
    ScaleMat2(w, ncol);
}

void Mat::ScaleMat(double* w, int m)
{
    Mat& S = (*this);
    for (register int i = 0; i < m; i++)
        if (row[i]) {
            for (register int j = 0; j < ncol; j++)
                if (col[j])
                    S(i, j) = S(i, j) * w[i];
        }
}

void Mat::ScaleMat2(double* w, int n)
{
    Mat& S = (*this);
    for (register int i = 0; i < nrow; i++)
        if (row[i])
        {
            for (register int j = 0; j < n; j++)
                if (col[j])
                    S(i, j) = S(i, j) * w[j];
        }
}


void Mat::ScaleMatInv(double* w)
{
    ScaleMatInv(w, nrow);
}

void Mat::ScaleMatInv2(double* w)
{
    ScaleMatInv2(w, ncol);
}

void Mat::ScaleMatInv(double* w, int m)
{
    Mat& S = (*this);
    for (register int i = 0; i < m; i++)
        if (row[i])
        {
            for (register int j = 0; j < ncol; j++)
                if (col[j])
                    S(i, j) = S(i, j) / w[i];
        }
}

void Mat::ScaleMatInv2(double* w, int n)
{
    Mat& S = (*this);
    for (register int i = 0; i < nrow; i++)
        if (row[i]) {
            for (register int j = 0; j < n; j++)
                if (col[j])
                    S(i, j) = S(i, j) / w[j];
        }
}

double Mat::getNorm()
{
    Mat& S = (*this);
    double ret = sqrt(S.Scal(S, S.nRow()));
    return ret;
};

void Mat::Print(char* s, bool m, int nfrac, int n)
{
    char format[24];
    Mat& S = *this;
    if (m)
        sprintf(format, "%s%d.%df ", "%", n, nfrac);
    else
        sprintf(format, "%s %s%d.%df ", "(%2d %2d)", "%", n, nfrac);

    if (m)
        printf("%s=[\n", s);
    else
        printf("%s: %s\n", cLable, s);
    for (register int i = 0; i < nrow; i++)
    {
        if (row[i])
        {
            if (m)
                printf("[");
            for (register int j = 0; j < ncol; j++)
                if (col[j])
                {
                    if (m)
                        printf(format, S(i, j));
                    else
                        printf(format, i, j, S(i, j));
                }
            if (m)
                printf("]\n");
            else
                printf("\n");
        }
    }
    if (m)
        printf("]\n");
}

void Mat::ShiftRow(int m)
{
    Mat& S = *this;
    for (register int j = 0; j < ncol; j++) {
        if (col[j]) {
            for (register int i = 0; i < nrow - m; i++) {
                S(i, j) = S(i + m, j);
                row[i] = row[i + m];
            }
        }
    }
    nrow -= m;
}

void Mat::SelfTransp()
{
    Mat& S = *this;
    for (register int i = 0; i < nrow; i++)
    {
        for (register int j = 0; j <= i; j++)
        {
            double tmp;
            tmp = S(j, i);
            S(j, i) = S(i, j);
            S(i, j) = tmp;
        }
    }

    for (register int i = 0; i < nrow; i++)
    {
        bool btmp = row[i];
        row[i] = col[i];
        col[i] = btmp;
    }
}

void Mat::Plus(const Mat& M1, const Mat& M2)
{
    Mat& S = *this;
    for (register int j = 0; j < M1.ncol; j++)
        col[j] = M1.col[j] && M2.col[j];
    for (register int i = 0; i < M1.nrow; i++)
    {
        row[i] = M1.row[i] && M2.row[i];
        if (row[i])
            for (register int j = 0; j < M1.ncol; j++)
                if (col[j])
                    S(i, j) = M1.val[i][j] + M2.val[i][j];
    }
    return;
};

void Mat::Minus(const Mat& M1, const Mat& M2)
{
    Mat& S = *this;
    for (register int j = 0; j < M1.ncol; j++)
        col[j] = M1.col[j] && M2.col[j];
    for (register int i = 0; i < M1.nrow; i++) {
        row[i] = M1.row[i] && M2.row[i];
        if (row[i])
            for (register int j = 0; j < M1.ncol; j++)
                if (col[j])
                    S(i, j) = M1.val[i][j] - M2.val[i][j];
    }
    return;
};

void Mat::Mult(const Mat& M1, const Mat& M2, double* w)
{
    Mult(false, 1., M1, M2, w);
};

void Mat::Mult_(const Mat& M1, const Mat& M2, double* w)
{
    Mult_(false, 1., M1, M2, w);
};

void Mat::TMult(const Mat& M1, const Mat& M2, double* w)
{
    TMult(false, 1., M1, M2, w);
};

void Mat::TMult(bool mode, double alp, const Mat& M1, const Mat& M2, double* w)
{
    Mat& S = *this;
    for (register int i = 0; i < M1.ncol; i++)
        row[i] = M1.col[i];
    for (register int j = 0; j < M2.ncol; j++)
        col[j] = M2.col[j];

    for (register int i = 0; i < M1.ncol; i++) {
        if (row[i]) {
            for (register int j = 0; j < M2.ncol; j++)
                if (col[j]) {
                    double sc = 0.;
                    for (register int k = 0; k < M1.nrow; k++)
                        if (M1.row[k] && M2.row[k])
                            sc += (M1.val[k][i] * M2.val[k][j] * (w ? w[k] : 1.));
                    if (mode)
                        S(i, j) += (alp * sc);
                    else
                        S(i, j) = (alp * sc);
                }
        }
    }
    return;
};

void Mat::Mult33(const Mat& M1, const Mat& M2)
{
    Mat& S = *this;
    S.val[0][0] = M1.val[0][0] * M2.val[0][0] + M1.val[0][1] * M2.val[1][0] + M1.val[0][2] * M2.val[2][0];
    S.val[0][1] = M1.val[0][0] * M2.val[0][1] + M1.val[0][1] * M2.val[1][1] + M1.val[0][2] * M2.val[2][1];
    S.val[0][2] = M1.val[0][0] * M2.val[0][2] + M1.val[0][1] * M2.val[1][2] + M1.val[0][2] * M2.val[2][2];

    S.val[1][0] = M1.val[1][0] * M2.val[0][0] + M1.val[1][1] * M2.val[1][0] + M1.val[1][2] * M2.val[2][0];
    S.val[1][1] = M1.val[1][0] * M2.val[0][1] + M1.val[1][1] * M2.val[1][1] + M1.val[1][2] * M2.val[2][1];
    S.val[1][2] = M1.val[1][0] * M2.val[0][2] + M1.val[1][1] * M2.val[1][2] + M1.val[1][2] * M2.val[2][2];

    S.val[2][0] = M1.val[2][0] * M2.val[0][0] + M1.val[2][1] * M2.val[1][0] + M1.val[2][2] * M2.val[2][0];
    S.val[2][1] = M1.val[2][0] * M2.val[0][1] + M1.val[2][1] * M2.val[1][1] + M1.val[2][2] * M2.val[2][1];
    S.val[2][2] = M1.val[2][0] * M2.val[0][2] + M1.val[2][1] * M2.val[1][2] + M1.val[2][2] * M2.val[2][2];
};

void Mat::TMult33(const Mat& M1, const Mat& M2)
{
    Mat& S = *this;
    S.val[0][0] = M1.val[0][0] * M2.val[0][0] + M1.val[1][0] * M2.val[1][0] + M1.val[2][0] * M2.val[2][0];
    S.val[0][1] = M1.val[0][0] * M2.val[0][1] + M1.val[1][0] * M2.val[1][1] + M1.val[2][0] * M2.val[2][1];
    S.val[0][2] = M1.val[0][0] * M2.val[0][2] + M1.val[1][0] * M2.val[1][2] + M1.val[2][0] * M2.val[2][2];

    S.val[1][0] = M1.val[0][1] * M2.val[0][0] + M1.val[1][1] * M2.val[1][0] + M1.val[2][1] * M2.val[2][0];
    S.val[1][1] = M1.val[0][1] * M2.val[0][1] + M1.val[1][1] * M2.val[1][1] + M1.val[2][1] * M2.val[2][1];
    S.val[1][2] = M1.val[0][1] * M2.val[0][2] + M1.val[1][1] * M2.val[1][2] + M1.val[2][1] * M2.val[2][2];

    S.val[2][0] = M1.val[0][2] * M2.val[0][0] + M1.val[1][2] * M2.val[1][0] + M1.val[2][2] * M2.val[2][0];
    S.val[2][1] = M1.val[0][2] * M2.val[0][1] + M1.val[1][2] * M2.val[1][1] + M1.val[2][2] * M2.val[2][1];
    S.val[2][2] = M1.val[0][2] * M2.val[0][2] + M1.val[1][2] * M2.val[1][2] + M1.val[2][2] * M2.val[2][2];
};

void Mat::Plus33(const Mat& M1, const Mat& M2)
{
    Mat& S = *this;
    S.val[0][0] = M1.val[0][0] + M2.val[0][0];
    S.val[0][1] = M1.val[0][1] + M2.val[0][1];
    S.val[0][2] = M1.val[0][2] + M2.val[0][2];

    S.val[1][0] = M1.val[1][0] + M2.val[1][0];
    S.val[1][1] = M1.val[1][1] + M2.val[1][1];
    S.val[1][2] = M1.val[1][2] + M2.val[1][2];

    S.val[2][0] = M1.val[2][0] + M2.val[2][0];
    S.val[2][1] = M1.val[2][1] + M2.val[2][1];
    S.val[2][2] = M1.val[2][2] + M2.val[2][2];
};

void Mat::Mult(bool mode, double alp, const Mat& M1, const Mat& M2, double* w)
{
    Mat& S = *this;
    bool lin = true;
    for (register int i = 0; i < M1.nrow; i++) {
        row[i] = M1.row[i];
        lin &= row[i];
    }
    for (register int j = 0; j < M2.ncol; j++) {
        col[j] = M2.col[j];
        lin &= col[j];
    }
    if (!w && (alp == 1.) && !mode) {
        double sc = 0.;
        if (lin) {
            for (register int i = 0; i < M1.nrow; i++)
                for (register int j = 0; j < M2.ncol; j++) {
                    sc = 0.;
                    for (register int k = 0; k < M1.ncol; k++)
                        sc += M1.val[i][k] * M2.val[k][j];
                    S(i, j) = sc;
                }
        }
        else {
            for (register int i = 0; i < M1.nrow; i++) {
                if (row[i]) {
                    for (register int j = 0; j < M2.ncol; j++)
                        if (col[j]) {
                            sc = 0.;
                            for (register int k = 0; k < M1.ncol; k++)
                                if (M1.col[k] && M2.row[k])
                                    sc += M1.val[i][k] * M2.val[k][j];
                            S(i, j) = sc;
                        }
                };
            };
        };
        return;
    }
    for (register int i = 0; i < M1.nrow; i++) {
        if (row[i]) {
            for (register int j = 0; j < M2.ncol; j++)
                if (col[j]) {
                    double sc = 0.;
                    for (register int k = 0; k < M1.ncol; k++)
                        if (M1.col[k] && M2.row[k])
                            sc += (M1.val[i][k] * M2.val[k][j] * (w ? w[k] : 1.));
                    if (mode)
                        S(i, j) += (alp * sc);
                    else
                        S(i, j) = (alp * sc);
                }
        }
    }
    return;
};

void Mat::Mult_(bool mode, double alp, const Mat& M1, const Mat& M2, double* w)
{
    Mat& S = *this;
    for (register int i = 0; i < M1.nrow; i++)
        row[i] = M1.row[i];
    for (register int j = 0; j < M2.ncol; j++)
        col[j] = M2.col[j];
    for (register int i = 0; i < M1.nrow; i++) {
        for (register int j = 0; j < M2.ncol; j++) {
            double sc = 0.;
            for (register int k = 0; k < M1.ncol; k++)
                sc += (M1.val[i][k] * M2.val[k][j] * (w ? w[k] : 1.));
            if (mode)
                S(i, j) += (alp * sc);
            else
                S(i, j) = (alp * sc);
        }
    }
    return;
};

void Mat::Mult(bool mode, double alp, LMat& M1, const Mat& M2)
{
    Mat& S = *this;
    int n = M1.nInd();

    for (register int i = 0; i < n; i++)
        row[i] = M1.Ind(i);

    for (register int j = 0; j < M2.ncol; j++)
        col[j] = M2.col[j];

    for (register int i = 0; i < n; i++) {
        if (row[i]) {
            for (register int j = 0; j < M2.ncol; j++)
                if (col[j]) {
                    double sc = 0.;
                    for (register int k = 0; k < n; k++)
                        if (M1.Ind(k) && M2.row[k]) {
                            double m1 = 0.;
                            if (i >= k)
                                m1 = M1(i, k);
                            else
                                m1 = M1(k, i);
                            sc += (m1 * M2.val[k][j]);
                        }
                    if (mode)
                        S(i, j) += (alp * sc);
                    else
                        S(i, j) = (alp * sc);
                }
        }
    }
    return;
};

/*
void Mat::operator=(const Mat & M)
{
  this->Move(M);
  //  printf(" =: %d\n", M.myNum_);
  ReleaseBuf(M.myNum_);
  //  printf(" =: isfree %d\n", isfree_[M.myNum_]);
}

void Mat::operator+=(const Mat &M)
{
  this->Plus(*this, M);
  //  printf(" +=: %d\n", M.myNum_);
  ReleaseBuf(M.myNum_);
  //  printf(" +=: isfree %d\n", isfree_[M.myNum_]);
}

Mat& operator+(const Mat &M1, const Mat &M2)
{
  //  Mat& T = *Mat::Buf();
  int ibuf = FindFree();
  //  printf(" + %d\n", ibuf);
  Mat& T = *Mat::Buf_[ibuf];
  T.Reset(M1.nrow, M1.ncol);
  T.Plus(M1, M2);
  return T;
};

void Mat::operator-=(const Mat &M)
{
  //  printf(" -=: %d\n", M.myNum_);
  ReleaseBuf(M.myNum_);
  //  printf(" -=: isfree %d\n", isfree_[M.myNum_]);
  this->Minus(*this, M);
}

Mat& operator-(const Mat &M1, const Mat &M2)
{
  int ibuf = FindFree();
  //  printf(" - %d\n", ibuf);
  Mat& T = *Mat::Buf_[ibuf];
  T.Reset(M1.nrow, M1.ncol);
  T.Minus(M1, M2);
  return T;
};

Mat& operator*(const Mat &M1, const Mat &M2)
{
  int ibuf = FindFree();
  //  printf(" * %d\n", ibuf);
  Mat& T = *Mat::Buf_[ibuf];
  T.Reset(M1.nrow, M2.ncol);
  T.Mult(M1, M2);
  return T;
};

Mat& operator&(const Mat &M1, const Mat &M2)
{
  int ibuf = FindFree();
  //  printf(" & %d\n", ibuf);
  Mat& T = *Mat::Buf_[ibuf];
  T.Reset(M1.ncol, M2.ncol);
  T.TMult(M1, M2);
  return T;
};
*/
void Mat::Transp(Mat& M)
{
    Mat& S = *this;
    for (register int j = 0; j < M.ncol; j++)
        row[j] = M.col[j];
    for (register int i = 0; i < M.nrow; i++)
    {
        col[i] = M.row[i];
        if (col[i])
            for (register int j = 0; j < M.ncol; j++)
                if (row[j])
                    S(j, i) = M(i, j);
    }
    return;
};

/*
Mat& operator~(const Mat &M)
{
  int ibuf = FindFree();
  //  printf(" ~ %d\n", ibuf);
  Mat& T = *Mat::Buf_[ibuf];
  T.Reset(M.ncol,M.nrow);
  T.Transp((Mat&)M);
  return T;
};
*/

#if MAT_DEBUG == 1
double& LMat::operator()(int i, int j)
{
    if (i < 0 || i >= nind) {
        printf(" LMAT out of range1 %d %d\n", i, nind);
        exit(999);
    }
    if (j < 0 || j > i) {
        printf(" LMAT out of range2 %d %d\n", i, j);
        exit(1000);
    }
    return val[i][j];
}

double*& LMat::operator [](int i)
{
    if (i < 0 || i >= nind) {
        printf(" LMAT out of range4 %d %d\n", i, nind);
        exit(1001);
    }
    return val[i];
}
#endif

LMat::LMat(int nind_, char cLable_[IdentLength])
{
    myType = HEAP;
    nind = maxind = nind_;
    dEpsilon = 1e-8;
    dEpsilonReg = 5e-8;
    strcpy(cLable, cLable_);
    storage = new double[nind * (nind + 1) / 2];
    ind = new bool[nind];
    for (register int k = 0; k < nind; k++)
        ind[k] = true;
    val = new double* [nind];
    for (register int k = 0; k < nind; k++)
        val[k] = storage + k * (k + 1) / 2;
    myNum_ = -1;
}

LMat::LMat(M_TYPE mt, int nind_, char cLable_[IdentLength])
{
    myType = mt;
    nind = maxind = nind_;
    dEpsilon = 1e-8;
    dEpsilonReg = 5e-8;
    strcpy(cLable, cLable_);
    if (mt == HEAP) {
        storage = new double[nind * (nind + 1) / 2];
        ind = new bool[nind];
    }
    if (mt == MALLOC) {
        storage = (double*)malloc(sizeof(double) * nind * (nind + 1) / 2);
        ind = (bool*)malloc(sizeof(bool) * nind);
    }
    for (register int k = 0; k < nind; k++)
        ind[k] = true;
    if (mt == HEAP)
        val = new double* [nind];
    if (mt == MALLOC)
        val = (double**)malloc(sizeof(double*) * nind);
    for (register int k = 0; k < nind; k++)
        val[k] = storage + k * (k + 1) / 2;
    myNum_ = -1;
}

LMat::LMat(const LMat& M)
{
    nind = M.nind;
    dEpsilon = M.dEpsilon;
    dEpsilonReg = M.dEpsilonReg;
    strcpy(cLable, M.cLable);
    storage = new double[nind * (nind + 1) / 2];
    memcpy((void*)storage, (void*)M.storage, sizeof(double) * nind * (nind + 1) / 2);
    ind = new bool[nind];
    for (register int k = 0; k < nind; k++)
        ind[k] = M.ind[k];

    val = new double* [nind];
    for (register int k = 0; k < nind; k++)
        val[k] = storage + k * (k + 1) / 2;
};

LMat::~LMat()
{
    //  printf("LM%d\n", myType);
    if (myType == HEAP) {
        delete[]val;
        delete[]storage;
        delete[]ind;
    }
    if (myType == MALLOC) {
        if (val)
            free(val);
        if (storage)
            free(storage);
        if (ind)
            free(ind);
    }
}

void LMat::Reset(int nind_)
{
    nind = nind_;
    for (register int k = 0; k < nind; k++)
    {
        ind[k] = true;
        val[k] = storage + k * (k + 1) / 2;
    }
}

void LMat::ResetAvail()
{
    for (register int k = 0; k < nind; k++) ind[k] = true;
}

void LMat::Move(const LMat& M)
{
    for (register int i = 0; i < M.nind; i++)
        ind[i] = M.ind[i];
    strcpy(cLable, M.cLable);
    memcpy((void*)storage, (void*)M.storage, sizeof(double) * M.nind * (M.nind + 1) / 2);
};

//void LMat::operator=(const LMat & M)
//{
//  this->Move(M);
//  ReleaseLBuf(M.myNum_);
//}

void LMat::Zero()
{
    LMat& S = *this;
    for (register int i = 0; i < nind; i++)
        if (ind[i])
        {
            for (register int j = 0; j <= i; j++)
                if (ind[j])
                    S(i, j) = 0.;
        }
}

void LMat::Ident()
{
    LMat& S = *this;
    S.Zero();
    for (register int i = 0; i < nind; i++)
        if (ind[i])
            S(i, i) = 1.;
}

void LMat::MultScal(double a)
{
    LMat& S = *this;
    for (register int i = 0; i < nind; i++)
    {
        if (ind[i])
        {
            for (register int j = 0; j <= i; j++)
            {
                if (ind[j])
                    S(i, j) = S(i, j) * a;
            }
        }
    }
}

void LMat::Print(char* s, bool m, int nfrac, int n)
{
    char format[24];
    LMat& S = *this;
    //  sprintf(format,"%s %s%d.%df ", "(%2d %2d)", "%", n, nfrac);
    //  printf("%s: %s \r\n", cLable, s);
    if (m)
        sprintf(format, "%s%d.%df ", "%", n, nfrac);
    else
        sprintf(format, "%s %s%d.%df ", "(%2d %2d)", "%", n, nfrac);

    if (m)
        printf("%s=[\n", s);
    else
        printf("%s: %s\n", cLable, s);

    for (register int i = 0; i < nind; i++)
        if (ind[i])
        {
            if (m)
                printf("[");
            for (register int j = 0; j < nind; j++)
                if (ind[j])
                {
                    if (m)
                        printf(format, ((i >= j) ? S(i, j) : S(j, i)));
                    else
                        printf(format, i, j, ((i >= j) ? S(i, j) : S(j, i)));
                }
            if (m)
                printf("]\n");
            else
                printf("\n");
        }
    if (m)
        printf("]\n");
}

void LMat::SymProd_AAt(double dAdd, Mat& A, double* w)
{
    LMat& S = *this;
    int N = A.nCol();
    int M = A.nRow();

    //  if(dAdd == 0.)
    //  {
    //    for(int i=0; i<M; i++)
    //      ind[i] = A.Row(i);
    //  }

    for (register int i = 0; i < M; i++) {
        if (ind[i] && A.Row(i)) {
            for (register int j = 0; j <= i; j++) {
                if (ind[j] && A.Row(j)) {
                    double sc = 0;
                    for (register int k = 0; k < N; k++) {
                        if (A.Col(k))
                            sc += A(i, k) * A(j, k) * (w ? w[k] : 1.);
                    }
                    S(i, j) += (dAdd * sc);
                }
            }
        }
    }
}

void LMat::SymProd_AtA(double dAdd, Mat& A, double* w)
{
    LMat& S = *this;
    int N = A.nCol();
    int M = A.nRow();
    //  if(dAdd == 0.)
    //  {
    //    for(int i=0; i<N; i++)
    //      ind[i] = A.Col(i);
    //  }

    for (register int i = 0; i < N; i++)
    {
        if (ind[i] && A.Col(i))
        {
            for (register int j = 0; j <= i; j++)
            {
                if (ind[j] && A.Col(j))
                {
                    double sc = 0;
                    for (register int k = 0; k < M; k++)
                    {
                        if (A.Row(k)) {
                            sc += A(k, i) * A(k, j) * (w ? w[k] : 1.);
                        }
                    }
                    //          printf("%d %d %f\n", i, j, sc);
                    S(i, j) += (dAdd * sc);
                }
            }
        }
    }
}

void LMat::SymProd_AAt(double dAdd, LMat& A, double* w)
{
    LMat& S = *this;
    int N = A.nInd();

    //  if(dAdd == 0.)
    //  {
    //    for(int i=0; i<N; i++)
    //      ind[i] = A.Col(i);
    //  }

    for (register int i = 0; i < N; i++) {
        if (ind[i] && A.Ind(i)) {
            for (register int j = 0; j <= i; j++) {
                if (ind[j] && A.Ind(j)) {
                    double sc = 0;
                    int ke = i;
                    if (j < i) ke = j;
                    for (register int k = 0; k <= ke; k++) {
                        if (A.Ind(k))
                            sc += A(i, k) * A(j, k) * (w ? w[k] : 1.);
                    }
                    S(i, j) += (dAdd * sc);
                }
            }
        }
    }
}

void LMat::SymProd_AtA(double dAdd, LMat& A, double* w)
{
    LMat& S = *this;
    int N = A.nInd();

    //  if(dAdd == 0.)
    //  {
    //    for(int i=0; i<N; i++)
    //      ind[i] = A.Col(i);
    //  }

    for (register int i = 0; i < N; i++)
    {
        if (ind[i] && A.Ind(i))
        {
            for (register int j = 0; j <= i; j++)
            {
                if (ind[j] && A.Ind(j))
                {
                    double sc = 0;
                    int kb = i;
                    if (j > i) kb = j;
                    for (register int k = kb; k < N; k++)
                    {
                        if (A.Ind(k))
                            sc += A(k, i) * A(k, j) * (w ? w[k] : 1.);
                    }
                    S(i, j) += (dAdd * sc);
                }
            }
        }
    }
}

void LMat::QRfactorization(Mat& Q, Mat& Abuf, int ireg, bool& ret)
{
    // A*Q = [L 0]
    int N = Abuf.nCol();
    int M = Abuf.nRow();
    nind = M;
    LMat& L = *this;

    for (register int i = 0; i < M; i++)
        ind[i] = true;

    Q.Reset(N, N);
    Q.Ident();
    for (register int i = 0; i < N; i++) {
        Q.Col(i) = Abuf.Col(i);
        Q.Row(i) = Abuf.Col(i);
    }

    double h[MAX_DIM];
    int ik = 0;
    int M_ = M;
    if (N < M_)
        M_ = N;
    for (register int k = 0; k < M; k++) {
        while (!Abuf.Col(ik))
            ik++;
        if (ik >= N)
            continue;
        for (register int i = 0; i < N; i++)
            h[i] = 0.;
        double sq = 0.;
        for (register int i = ik; i < N; i++) {
            if (!Abuf.Col(i))
                continue;
            h[i] = Abuf(k, i);
            sq += (h[i] * h[i]);
        }
        double alp = -sqrt(sq);
        double sq_ha = sq + h[ik] * alp;
        if (fabs(sq_ha) < dEpsilon) {
            ik++;
            continue;
        }
        double tau = 0.5 / (sq_ha);
        h[ik] += alp;
        double Qh[MAX_DIM];
        for (register int i = 0; i < N; i++)
            Qh[i] = 0.;
        for (register int i = 0; i < N; i++)
        {
            if (!Abuf.Col(i))
                continue;
            Qh[i] = 0.;
            for (register int j = 0; j < N; j++)
            {
                if (!Abuf.Col(j))
                    continue;
                Qh[i] += Q(i, j) * h[j];
            }
        }
        for (register int i = 0; i < N; i++)
        {
            if (!Abuf.Col(i))
                continue;
            for (register int j = 0; j < N; j++)
            {
                if (!Abuf.Col(j))
                    continue;
                Q(i, j) -= (2. * tau * Qh[i] * h[j]);
            }
        }

        for (register int l = k; l < M; l++)
        {
            double sqa = 0.;
            for (register int j = 0; j < N; j++)
            {
                if (!Abuf.Col(j))
                    continue;
                sqa += (Abuf(l, j) * h[j]);
            }
            for (register int j = 0; j < N; j++)
            {
                if (!Abuf.Col(j))
                    continue;
                Abuf(l, j) -= (2. * tau * sqa * h[j]);
            }
            L(l, k) = Abuf(l, ik);
        }
        //    Abuf.Print(" Abuf ", 4);
        ik++;
    }
    ik = 0;
    ret = true;
    for (register int k = 0; k < M; k++)
    {
        while (!Abuf.Col(ik))
            ik++;
        if (ik >= N)
            continue;
        for (register int l = k; l < M; l++)
            L(l, k) = Abuf(l, ik);
        if ((k >= ireg) && (fabs(L(k, k)) <= dEpsilon))
            L(k, k) = 2. * dEpsilonReg;
        ret = ret && (fabs(L(k, k)) > dEpsilon);
        //    if(!ret)
        //      printf("k %d L %12f\n", k, L(k,k));
        ik++;
    }
};

void LMat::Chol(LMat& L, bool& ret, int iStart, int iStop, int mark)
{
    (void)mark;
    double tmp, sq;
    double* Li, * Lj;
    double* Lik, * Ljk;
    int ir, jr;
    bool* indc;

    ret = true;
    for (register int i = iStart; i <= iStop; i++)
    {
        L.ind[i] = ind[i];
        Li = L.val[i]; Lj = val[i];
        for (register int j = iStart; j <= i; j++, Li++, Lj++)
            *Li = *Lj;
    }
    ir = 0;
    indc = ind + iStart;
    for (register int i = iStart; i <= iStop; i++)
    {
        if (*indc)
        {
            Li = L.val[i];
            tmp = *(Li + i);
            //      printf("i %d tmp %f \n", i, tmp); 
            if (ir > 0)
            {
                jr = 0;
                for (register int j = iStart; j < i; j++)
                {
                    if (ind[j])
                    {
                        sq = 0.;
                        Lj = L.val[j];
                        if (jr > 0)
                        {
                            Lik = Li; Ljk = Lj;
                            for (register int k = iStart; k < j; k++)
                            {
                                if (ind[k])
                                {
                                    sq = sq + *Lik * *Ljk;
                                };
                                Lik++; Ljk++;
                            }
                        }
                        jr++;
                        if (fabs(*(Lj + j)) < dEpsilon)
                        {
                            //              printf("!Chol %f j %d\r\n", *(Lj+j), j); 
                            ret = false; return;
                        }
                        sq = (*(Li + j) - sq);
                        *(Li + j) = sq / *(Lj + j);
                        tmp = tmp - *(Li + j) * *(Li + j);
                        //	    printf("j %d tmp %f\r\n", j, tmp); 
                    };
                };
            };
            ir++;
            if (tmp < dEpsilon) {
                //        printf("!!Chol %f i %d\r\n", tmp, i); 
                ret = false;
                return;
            };
            *(Li + i) = sqrt(tmp);
        };
        indc++;
    }
}

void LMat::Chol_(LMat& L, bool& ret, int iStart, int iStop, int mark)
{
    double tmp, sq;
    double* Li, * Lj;
    double* Lik, * Ljk;
    int ir, jr;
    bool* indc;

    ret = true;
    for (register int i = iStart; i <= iStop; i++)
    {
        L.ind[i] = ind[i];
        Li = L.val[i]; Lj = val[i];
        for (register int j = iStart; j <= i; j++, Li++, Lj++)
            *Li = *Lj;
    }
    ir = 0;
    indc = ind + iStart;
    for (register int i = iStart; i <= iStop; i++) {
        Li = L.val[i];
        tmp = *(Li + i);
        //      printf("i %d tmp %f \n", i, tmp); 
        if (ir > 0) {
            jr = 0;
            for (register int j = iStart; j < i; j++) {
                sq = 0.;
                Lj = L.val[j];
                if (jr > 0)
                {
                    Lik = Li; Ljk = Lj;
                    for (register int k = iStart; k < j; k++) {
                        {
                            sq = sq + *Lik * *Ljk;
                        };
                        Lik++; Ljk++;
                    }
                }
                jr++;
                if (fabs(*(Lj + j)) < dEpsilon) {
                    printf("!Chol %f j %d mark %d\r\n", *(Lj + j), j, mark);
                    ret = false; return;
                }
                sq = (*(Li + j) - sq);
                *(Li + j) = sq / *(Lj + j);
                tmp = tmp - *(Li + j) * *(Li + j);
                //	    printf("j %d tmp %f\r\n", j, tmp); 
            };
        };
        ir++;
        if (tmp < dEpsilon) {
            printf("!!Chol %f mark %d i %d\r\n", tmp, mark, i);
            ret = false;
            return;
        };
        *(Li + i) = sqrt(tmp);
    };
    indc++;
}

void LMat::Forward(Mat& M, int iStart)
{
    LMat& L = (*this);
    double sq;
    int i, j, k;
    for (k = 0; k < M.nCol(); k++)
    {
        if (M.Col()[k])
        {
            for (i = 0; i < (nind - iStart); i++)
            {
                if (ind[i])
                {
                    sq = 0.;
                    for (j = 0; j < i; j++)
                    {
                        if (ind[j])
                            sq = sq + L(iStart + i, iStart + j) * M(j, k);
                    }
                    M(i, k) = M(i, k) - sq;
                    M(i, k) = M(i, k) / L(iStart + i, iStart + i);
                }
            }
        }
    }
}

void LMat::Forward_(Mat& M, int iStart)
{
    LMat& L = (*this);
    double sq;
    int i, j, k;
    for (k = 0; k < M.nCol(); k++) {
        for (i = 0; i < (nind - iStart); i++) {
            sq = 0.;
            for (j = 0; j < i; j++)
                sq = sq + L(iStart + i, iStart + j) * M(j, k);
            M(i, k) = M(i, k) - sq;
            M(i, k) = M(i, k) / L(iStart + i, iStart + i);
        }
    }
}

void LMat::Backward(Mat& M)
{
    LMat& L = (*this);
    for (register int l = 0; l < M.nCol(); l++)
        if (M.Col(l))
            for (register int k = 0; k < nind; k++)
            {
                int i = nind - k - 1;
                if (ind[i])
                {
                    M(i, l) /= L(i, i);
                    for (register int j = 0; j < i; j++)
                        if (ind[j])
                            M(j, l) -= L(i, j) * M(i, l);
                }
            }
    return;
}

void LMat::Backward_(Mat& M)
{
    LMat& L = (*this);
    for (register int l = 0; l < M.nCol(); l++) {
        for (register int k = 0; k < nind; k++) {
            int i = nind - k - 1;
            M(i, l) /= L(i, i);
            for (register int j = 0; j < i; j++)
                M(j, l) -= L(i, j) * M(i, l);
        }
    }
    return;
}

void LMat::Plus(const LMat& M1, const LMat& M2)
{
    LMat& L = (*this);
    for (register int j = 0; j < M1.nind; j++)
        ind[j] = M1.ind[j] && M2.ind[j];
    for (register int i = 0; i < M1.nind; i++)
        if (ind[i])
            for (register int j = 0; j <= i; j++)
                if (ind[j])
                    L(i, j) = M1.val[i][j] + M2.val[i][j];
    return;
};

void LMat::operator+=(const LMat& M)
{
    this->Plus(*this, M);
}

/*
LMat& operator+(const LMat &M1, const LMat &M2)
{
  int ibuf = FindLFree();
  LMat& T = *LMat::Buf_[ibuf];
  T.Reset(M1.nind);
  T.Plus(M1, M2);
  return T;
  //  LMat& T = LMat::Buf();
  //  T.Reset(M1.nind);
  //  T.Plus(M1, M2);
  //  return T;
};
*/
void LMat::Lswap(int ia)
{
    LMat& T = *this;
    double dn2, dn, sc, tau, vii;
    double sw, Lii_;
    static double* Ljj_;
    //  int is;

    for (register int is = 0; is <= ia; is++) {
        sw = T(ia, is);
        T(ia, is) = T(ia + 1, is);
        T(ia + 1, is) = sw;
    };
    Lii_ = T(ia + 1, ia + 1);
    T(ia + 1, ia + 1) = 0.;

    dn2 = sq(T(ia, ia)) + sq(Lii_);
    dn = sqrt(dn2);
    tau = dn2 + T(ia, ia) * dn;
    vii = T(ia, ia) + dn;
    for (register int k = ia + 1; k < nind; k++) {
        sc = (T(ia, ia) * T(k, ia) + Lii_ * T(k, ia + 1)
            + dn * T(k, ia)) / tau;
        T(k, ia) = -(T(k, ia) - sc * vii);
        T(k, ia + 1) = T(k, ia + 1) - sc * Lii_;
        if (k > (ia + 1)) {
            Ljj_ = this->val[ia + 1] + ia + 1;
            if (*Ljj_ < 0.)
                T(k, ia + 1) = -T(k, ia + 1);
        }
    }
    T(ia, ia) = dn;
    if (T(ia + 1, ia + 1) < 0.)
        T(ia + 1, ia + 1) = -T(ia + 1, ia + 1);
}

#define ALP_MOD 0.01

bool LMat::ModLFact(Mat& a, int col, double alpha, Mat& p, Mat& w, bool bpr, bool breg)
{
    double e2 = 1.e-12;//dEpsilon * dEpsilon;
    int n = nInd();
    if (col < 0 || col >= a.nCol())
    {
        if (bpr)
            printf("!Modl col %d\r\n", col);
        return false;
    }
    int i;

    p.Reset(n, 1);
    w.Reset(n, 1);
    LMat& S = (*this);
    if (bpr)
        printf("Modl n %d col %d alpha %f\r\n", n, col, alpha);

    for (i = 0; i < n; i++)
    {
        if (Ind(i))
            p(i) = a(i, col);
    }
    S.Forward(p);
    double s = 0.;
    for (i = 0; i < n; i++)
    {
        if (Ind(i))
        {
            w(i) = a(i, col);
            s += sq(p(i));
        }
    }
    double D = 1. + alpha * s;
    if (D < e2) {
        //    printf("modl D %.18f\n", D);
        if (breg) {
            if (D > -e2) {
                D = ALP_MOD;
                //	printf("modl corr D %.18f\n", D);
                return false;
            }
            else
                return false;
        }
        else
            return false;
    }

    double sigma = alpha / (1. + sqrt(D));

    for (i = 0; i < n; i++)
    {
        if (Ind(i))
        {
            double q = sq(p(i));
            double theta = 1. + sigma * q;
            s -= q;
            double rho2 = sq(theta) + sq(sigma) * q * s;
            double rho = sqrt(rho2);
            S(i, i) = S(i, i) * rho;
            alpha /= rho2;
            sigma *= (1. + rho) / (rho * theta + rho2);
            for (register int k = i + 1; k < n; k++)
            {
                if (Ind(k))
                {
                    w(k) = w(k) - p(i) * S(k, i);
                    (*this)(k, i) = rho * (S(k, i) + alpha * p(i) * w(k));
                }
            }
        }
    }
    return true;
}

void LMat::SetSub(int ib, int jb, Mat& M,
    int imb, int ime, int jmb, int jme)
{
    LMat& S = *this;
    for (register int i = imb; i <= ime; i++)
        for (register int j = jmb; j <= jme; j++)
            S(ib + i - imb, jb + j - jmb) = M(i, j);
};

void LMat::SetSub(int ib, int jb, LMat& M,
    int imb, int ime, int jmb, int jme)
{
    LMat& S = *this;
    for (register int i = imb; i <= ime; i++)
        for (register int j = jmb; j <= jme; j++)
            S(ib + i - imb, jb + j - jmb) = M(i, j);
};

void LMat::SetSub(int ib, LMat& M,
    int imb, int ime)
{
    LMat& S = *this;
    for (register int i = imb; i <= ime; i++)
        for (register int j = imb; j <= i; j++)
            S(ib + i - imb, ib + j - imb) = M(i, j);
};
