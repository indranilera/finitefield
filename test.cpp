#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<string>
#include "finiteField.h"
using namespace std;


int kernel_main()
{
  fflush(stdout);
  Field* f = new Field(2,6,(unsigned char*) "101001");// ir = 1+x^2+x^5
  //
  int r,c;
  r = 2;c=5;
  unsigned char** a = new unsigned char*[r];
  for(int i = 0;i<r;i++)
    {
      a[i] = new unsigned char[c];
      for(int j = 0; j < c;j++)
	{
	  a[i][j]= 1;
	}
    }
  
  //for 4
  //a[1][0] = 0;a[1][1] = 0;a[1][2] = 0;
  //a[2][0] = 0;a[2][1] = 0;a[2][2] = 0;
  
  // for 5
  a[0][1] = 0;a[0][4] = 0;
  a[1][0] = 0;a[1][3]=  0;
  /////////////////////////////////////////
  for(int i = 0; i< r;i++)
    {
    for(int j = 0;j<c;j++)
      {
	int vv = (int)a[i][j];
	cout<<vv<<"   ";
      }
    cout<<endl;
    }
  //cout<<"\n--------------------\n";
  //////////////////////////////////////////
  Element** basis = 0;
  int rows = 0;
  f->GetKarnel(a,r,c,&basis,rows);
  //
  cout<<"\nthe basis : \n\n";
  for(int i = 0 ; i < rows;i++)
    {
      basis[i]->Show();cout<<endl;
    }
  cout<<endl;
  //
  for(int i = 0; i < rows;i++)
    {
      delete basis[i];
      basis[i] = 0;
    }
  delete [] basis;
  basis =0;
  //
  for(int i = 0; i < r ; i++)
    {
      delete [] a[i];
      a[i] = 0;
    }
  delete [] a;
  a = 0;
  //
  delete f;
  return 1;
}



int trace_main()
{
  fflush(stdout);
  Field* f = new Field(2,4,(unsigned char*) "1101");// ir = 1+x+x^3
  Element* a = new Element((unsigned char*)"010",3);
  //
  Element* b = new Element((unsigned char*)"100",3);
  //
  for(int i = 0;i<7;i++)
    {
      Element* c = 0;
      Element* d = 0;
      f->Power(a,i,&d);
      Element*e = 0;
      f->Multiply(d,b,&e);
      f->Trace(e,&c);
      c->Show();cout<<endl;
      //
      delete c;c=0;
      delete d;d=0;
      delete e;e = 0;
    }
  cout<<"\n\n\n";
  delete a;a=0;
  delete b;b=0;
  //delete c;c=0;
  delete f;f=0;
  return 0;
}


int SmallerCosetForC(Field* f,int kPower,int aPower)
{
  fflush(stdout);
  int fl = f->_len-1;
  //Field* f = new Field(2,fl+1,(unsigned char*) "1101");// ir = 1+x+x^3
  //Element* a = new Element((unsigned char*)"010",fl);
  //
  //Forming the matrix.....
  int r = fl;
  // Element*** mtrx  =0;// new Element**[r];
 //  f->MatrixAllocate(&mtrx,r,r);
 // //
 //  long int p =1;
 //  for(int i = 0; i< r;i++)
 //    {
 //      Element*c = 0;
 //      f->Power(f->_a,p,&c);
 //      p *= 2;
 //      for(int j = 0 ; j < r;j++)
 // 	{
 // 	  f->Power(c,j,&mtrx[j][i]);
 // 	}
 //      delete c;
 //      c=0;
 //    }
 
 //  // print the matrix
 //  //f->MatrixPrint(mtrx,r,r);
 //  //////  invert the matrix////////////////////
 //  Element*** invMatrix = 0;
 //  f->MatrixInvert(mtrx,r,&invMatrix);
 //  //
 //  Element*** chk = 0;
 //  int rr;
 //  int cc;
 //  //r = 3;
 //  //cout<<"i am here";
 //  f->MatrixMultiply(mtrx,r,r,invMatrix,r,r,&chk,rr,cc);
 //  //f->MatrixPrint(chk,rr,cc);
 //  /////////////////////////////////////////////
 //  Element*** dualBasis=0;
 //  f->MatrixAllocate(&dualBasis,r,1);
 //  //
 //  for(int i = 0; i<r;i++)
 //    {
 //      dualBasis[i][0] = new Element(invMatrix[0][i]);
 //    }
 //  //
 //  unsigned char** bMatrix = 0;
 //  f->ConvertToLowerFieldMatrix(dualBasis,r,1,&bMatrix);
 //  //
 //  unsigned char** bInv = 0;
 //  f->BinaryMatrixInvert(bMatrix,r,&bInv);
 //  //
 //  //cout<<"\n\nbinary matrix\n\n";
 //  for(int i = 0; i< r; i++)
 //    {
 //      for(int j = 0;j< r ;j++)
 // 	{
 // 	  int bv = bInv[i][j];
 // 	  //cout<<bv<<"    ";
 // 	}
 //      //cout<<endl;
 //    }

  Element*** traceCoff = 0;//new Element*[fl];
  //
  //cout<<"\nkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk\n";
  f->MatrixAllocate(&traceCoff,1,fl);
  long long int s = aPower;
  int m = f->GetCosetSize(s);
  for(int i = 0; i<fl;i++)
    {
      //
      Element* ap = 0;
      f->Power(f->_a,i,&ap);
      f->Trace(ap,&traceCoff[0][i],m);
      delete ap;
      ap=0;
    }
  //
  //cout<<"\ngggggggggggggggggggggggggg\n";
  unsigned char** binTraceMatrix = 0;
  f->ConvertToLowerFieldMatrix(traceCoff,1,fl,&binTraceMatrix);
  //
  Element** basis = 0;
  int rows = 0;
  f->GetKarnel(binTraceMatrix,fl,fl,&basis,rows);
  //
  cout<<"\nthe basis : \n\n";
  for(int i = 0 ; i < rows;i++)
    {
      basis[i]->Show();cout<<endl;
    }
  cout<<endl;

  // combine to get b0 
  // TO DO:
  for(int allCombination = 0; allCombination < rows;allCombination++)
    {
      noOfEqns++;
      Element * cValue = new Element(f->_zero);
      for(int i = 0; i<fl;i++)
	{
	  if(basis[allCombination]->_val[i] ==1)
	    {
	      Element* ap = 0;
	      f->Power(f->_a,i,&ap);
	      Element *bp = 0;
	      f->Add(ap,cValue,&bp);
	      delete cValue;
	      cValue = 0;
	      cValue = new Element(bp);
	      delete ap;ap=0;
	      delete bp;bp = 0;
	    }
	}
      ///////////////////////////////////////////
      // get vj's thereof
      unsigned char** cBinaryMatrix = 0;
      f->ConvertToLowerFieldMatrix(cValue,true,&cBinaryMatrix);
  
      unsigned char** vValue = 0;
      f->BinaryMatrixMultiply(cBinaryMatrix,1,r,f->_bInv,r,r,&vValue);
      //
      cout<<"\n\n v cofficient matrix Matrix vj \n\n";
      for(int i = 0 ; i < r ; i++)
	{  
	  int vvv = (int)vValue[0][i];
	  cout<<vvv<<"     ";
	}
      cout<<endl<<endl;
      //////////////////////////////////CHECK/////////////////////////////////////
      unsigned char bit = 0;
      for(int i = 0; i< -1;i++)//pow(2,fl)-1;i++)
	{
	  Element* ap = 0;
	  Element* apy = 0;
	  f->Power(f->_a,i,&ap);
	  f->Power(ap,aPower,&apy);
	  Element* traceOf =0;
	  f->Multiply(cValue,apy,&traceOf);
	  Element* traceValue = 0;
	  f->Trace(traceOf,&traceValue);
	  bit = traceValue->_val[0];
	  //
	  unsigned char s2 = 0;
	  for( int yj = 0; yj < r;yj++)
	    {
	      s2 =  s2 & ~(vValue[0][yj] & apy->_val[yj]) |
		(~s2 & (vValue[0][yj] & apy->_val[yj]));
	      s2 &= 1;
	    }
	  int v1 = (int)bit;
	  int v2 = (int)s2;
	  if(v1 != v2)
	    {
	      cout<<"\n\n its not matching \n\n";
	    }    
	  //
	  delete ap,apy,traceOf,traceValue;
	  ap=0;apy=0;traceOf = 0;traceValue = 0;
	}

      ////////////delete section//////////////////////////////
      delete cValue ;cValue = 0;
      
      //for(int i = 0;i<r;i++)
      //{
	  delete[] cBinaryMatrix[0];
	  cBinaryMatrix[0] = 0;
	  //}
      delete [] cBinaryMatrix;
      cBinaryMatrix = 0;
      //
      delete [] vValue[0];
      vValue[0] = 0;
      delete [] vValue;
      vValue = 0;
      
    }
  //**************************************************************//
for(int i = 0; i < rows;i++)
    {
      delete basis[i];
      basis[i] = 0;
    }
  delete [] basis;
  basis =0;

  //////////////////////////////////////////////////
  // delete the matrix
  //  f->MatrixDelete(&mtrx,r,r);
  // f->MatrixDelete(&invMatrix,r,r);
  //f->MatrixDelete(&chk,rr,cc);
  //f->MatrixDelete(&dualBasis,r,1);
  f->MatrixDelete(&traceCoff,1,fl);
  //
  // for(int i = 0;i<r;i++)
  //   {
  //     delete [] bMatrix[i];
  //     bMatrix[i] = 0;
  //     //
  //     delete [] bInv[i];
  //     bInv[i] = 0;
  //     //
  //   }
    
  delete []binTraceMatrix[0];
  binTraceMatrix[0] = 0;
  delete [] binTraceMatrix;
  binTraceMatrix = 0;
  //
  // delete [] bMatrix;
  // bMatrix = 0;
  // delete [] bInv;
  // bInv = 0;
  return 1;
}



int SmallerCosetForBk(Field*f,int kPower,int aPower)
{
  fflush(stdout);
   int fl = f->_len-1;
  //Field* f = new Field(2,fl+1,(unsigned char*) "1101");// ir = 1+x+x^3
  //Element* a = new Element((unsigned char*)"010",fl);
  //
  //Forming the matrix.....
  int r = fl;
  // Element*** mtrx  =0;// new Element**[r];
 //  f->MatrixAllocate(&mtrx,r,r);
 // //
 //  long int p =1;
 //  for(int i = 0; i< r;i++)
 //    {
 //      Element*c = 0;
 //      f->Power(f->_a,p,&c);
 //      p *= 2;
 //      for(int j = 0 ; j < r;j++)
 // 	{
 // 	  f->Power(c,j,&mtrx[j][i]);
 // 	}
 //      delete c;
 //      c=0;
 //    }
 
 //  // print the matrix
 //  //f->MatrixPrint(mtrx,r,r);
  
 //  //////  invert the matrix////////////////////
 //  Element*** invMatrix = 0;
 //  f->MatrixInvert(mtrx,r,&invMatrix);
 //  //cout<<"i am here";
 //  //f->MatrixPrint(invMatrix,r,r);
 //  //cout<<"i am here";

 //  //
 //  Element*** chk = 0;
 //  int rr;
 //  int cc;
 //  //r = 3;
 //  //cout<<"i am here";
 //  f->MatrixMultiply(mtrx,r,r,invMatrix,r,r,&chk,rr,cc);
 //  //f->MatrixPrint(chk,rr,cc);
 //  /////////////////////////////////////////////
 //  Element*** dualBasis=0;
 //  f->MatrixAllocate(&dualBasis,r,1);
 //  //
 //  for(int i = 0; i<r;i++)
 //    {
 //      dualBasis[i][0] = new Element(invMatrix[0][i]);
 //    }
 //  //
 //  unsigned char** bMatrix = 0;
 //  f->ConvertToLowerFieldMatrix(dualBasis,r,1,&bMatrix);
 //  //
 //  unsigned char** bInv = 0;
 //  f->BinaryMatrixInvert(bMatrix,r,&bInv);
 //  //
 //  cout<<"\n\nbinary matrix\n\n";
 //  for(int i = 0; i< r; i++)
 //    {
 //      for(int j = 0;j< r ;j++)
 // 	{
 // 	  int bv = bInv[i][j];
 // 	  cout<<bv<<"    ";
 // 	}
 //      cout<<endl;
 //    }
  /////////////////////////////////////////////
  // b0 b1 
  Element*** traceCoff = 0;//new Element*[fl];
  //
  //cout<<"\nkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk\n";
  f->MatrixAllocate(&traceCoff,1,fl);
  int s = pow(2,kPower) + aPower % (int) (pow(2,fl)-1);
  int m = f->GetCosetSize(s);
  cout<<"\n m = "<<m<<endl;
  for(int i = 0; i<fl;i++)
    {
      //
      Element* ap = 0;
      f->Power(f->_a,i,&ap);
      f->Trace(ap,&traceCoff[0][i],m);
      delete ap;
      ap=0;
    }
  //////////////////////////////////////
  

  /////////////////////////////////////
  //
  //cout<<"\ngggggggggggggggggggggggggg\n";
  unsigned char** binTraceMatrix = 0;
  f->ConvertToLowerFieldMatrix(traceCoff,1,fl,&binTraceMatrix);
  // cout<<"\ntrace matrix\n";
  // for(int i =0;i<fl;i++)
  //   {
  //     for(int j = 0 ; j < fl;j++)
  // 	{
  // 	  int vv = (int)binTraceMatrix[i][j];
  // 	  cout<<vv<<"   ";
  // 	}
  //     cout<<endl;
  //   }
  // cout<<endl;
  //////////////////////////////////////////
  //
  Element** basis = 0;
  int rows = 0;
  f->GetKarnel(binTraceMatrix,fl,fl,&basis,rows);
  //
  cout<<"\nthe basis : \n\n";
  for(int i = 0 ; i < rows;i++)
    {
      basis[i]->Show();cout<<endl;
    }
  cout<<endl;
  cout<<"\nrows = "<<rows<<endl;
  // combine to get b0 
  // TO DO:
  for(int allCombination = 0; allCombination < rows;allCombination++)
    {
      noOfEqns++;
      Element * b0 = new Element(f->_zero);
      for(int i = 0; i<fl;i++)
	{
	  if(basis[allCombination]->_val[i] ==1)
	    {
	      Element* ap = 0;
	      f->Power(f->_a,i,&ap);
	      Element *bp = 0;
	      f->Add(ap,b0,&bp);
	      delete b0;
	      b0 = 0;
	      b0 = new Element(bp);
	      delete ap;ap=0;
	      delete bp;bp = 0;
	    }
	}
      ////////////////////////////////////////////
      Element***bCol = 0;//new Element**[3];
      f->MatrixAllocate(&bCol,r,1);
      //
      
      //then form bCol
      for(int bi = 0 ; bi < r;bi++)
	{
	  if(bi == kPower)
	    {
	      bCol[bi][0] = new Element(b0);
	    }
	  else
	    {
	  bCol[bi][0] = new Element(f->_zero);
	    }
	}
      
      // after b0,b1 etc formed
      Element *** tMatrix = 0;
      int tr = 0;
      int tc = 0;
      f->MatrixMultiply(f->_mtrx,r,r,bCol,r,1,&tMatrix,tr,tc);
      unsigned char** tBinaryMatrix = 0;
      f->ConvertToLowerFieldMatrix(tMatrix,tr,tc,&tBinaryMatrix);
      //
      // cout<<"\n\n t binary cofficient matrix\n\n";
      // for(int i = 0; i< r; i++)
      //   {
      //     for(int j = 0;j< r;j++)
      // 	{
      // 	  int bv = tBinaryMatrix[i][j];
      // 	  cout<<bv<<"    ";
      // 	}
      //     cout<<endl;
      //   }
      
      ///////////////////////////////
      // to multiply tI
      unsigned char** aMatrix = 0;
      f->BinaryMatrixMultiply(tBinaryMatrix,r,r,f->_bInv,r,r,&aMatrix);
      cout<<"\n\nbinary cofficient matrix aij \n\n";
      for(int i = 0; i< r; i++)
	{
	  for(int j = 0;j< r;j++)
	    {
	      int bv = aMatrix[i][j];
	      cout<<bv<<"  ,  ";
	    }
	  cout<<endl;
	}
      
      /////////////////////CHECK/////////////////////////////////////
      unsigned char sum = 0;
      unsigned char xBit = 0;
      for(int p = 0 ; p < -1;p++)//pow(2,fl)-1;p++)
	{
	  //cout<<"p = "<<p<<endl;
	  //cout<<"\n-----------gggggggggggg---------------------------\n";
	  //a->Show();
	  Element* ap = 0;
	  f->Power(f->_a,p,&ap);
	  Element* apy=0;
	  if(aPower == -1)
	    {
	      f->Invert(ap,&apy);
	    }
	  else
	    {
	      //
	      f->Power(ap,aPower,&apy);
	    }
	  //ap->Show();
	  //xBit = ap->_val[bitToBeMatched];
	  s = 0;
	  for(int i = 0; i< r;i++)
	    {
	      for(int j = 0; j<r;j++)
		{
		  sum =  (sum & ~( aMatrix[i][j] & ap->_val[i] & apy->_val[j])) | (~sum &(aMatrix[i][j] & ap->_val[i] & apy->_val[j]));
		}
	    }
	  
	  int v1 = (int)xBit;
	  int v2 = (int)sum;
	  if(v1 == v2)
	    {
	      //cout<<"\n\n MATCHING.................\n\n";
	    }
	  else{
	    ap->Show();cout<<endl;apy->Show();
	    cout<<"\n\n not mqatching for p ="<<p<<"\n\n";
	  }
	  delete ap;
	  ap = 0;
	  delete apy;
	  apy=0;
	}
      //////////////////////////////// DELETE///////////////////////
      for(int i = 0;i<r;i++)
    {
      delete [] aMatrix[i];
      aMatrix[i] = 0;
      //
      //delete [] bMatrix[i];
      //bMatrix[i] = 0;
      //
      delete [] tBinaryMatrix[i];
      tBinaryMatrix[i] = 0;
      //
      //delete [] cBinaryMatrix[i];
      //cBinaryMatrix[i] = 0;
     
      //
      //delete [] bInv[i];
      //bInv[i] = 0;
    }
    
      // delete []binTraceMatrix[0];
      // binTraceMatrix[0] = 0;
      // delete [] binTraceMatrix;
      // binTraceMatrix = 0;
      delete [] aMatrix;
      aMatrix = 0;
      delete [] tBinaryMatrix;
      tBinaryMatrix = 0;
      f->MatrixDelete(&tMatrix,tr,tc);
      f->MatrixDelete(&bCol,r,1);
      delete b0;b0 = 0;    
    }

  //*********************************************************//
  for(int i = 0; i < rows;i++)
    {
      delete basis[i];
      basis[i] = 0;
    }
  delete [] basis;
  basis =0;

  //////////////////////////////////////////////////
  // delete the matrix
  // f->MatrixDelete(&mtrx,r,r);
  // f->MatrixDelete(&invMatrix,r,r);
  // f->MatrixDelete(&chk,rr,cc);
  // f->MatrixDelete(&dualBasis,r,1);
  f->MatrixDelete(&traceCoff,1,fl);
  //
  // for(int i = 0;i<r;i++)
  //   {
  //     delete [] bMatrix[i];
  //     bMatrix[i] = 0;
  //     //
  //     delete [] bInv[i];
  //     bInv[i] = 0;
  //     //
  //   }
    
  delete []binTraceMatrix[0];
  binTraceMatrix[0] = 0;
  delete [] binTraceMatrix;
  binTraceMatrix = 0;
  //
  // delete [] bMatrix;
  // bMatrix = 0;
  // delete [] bInv;
  // bInv = 0;
  return 1;
}


int SameBigCosetForBk(Field* f,int k1Power,int k2Power,int aPower)
{
  fflush(stdout);
  int fl = f->_len-1;
  //Field* f = new Field(2,fl+1,(unsigned char*) "1101");// ir = 1+x+x^3
  //Element* a = new Element((unsigned char*)"010",fl);
  //
  //Forming the matrix.....
  int r = fl;
  // Element*** mtrx  =0;// new Element**[r];
 //  f->MatrixAllocate(&mtrx,r,r);
 // //
 //  long int p =1;
 //  for(int i = 0; i< r;i++)
 //    {
 //      Element*c = 0;
 //      f->Power(f->_a,p,&c);
 //      p *= 2;
 //      for(int j = 0 ; j < r;j++)
 // 	{
 // 	  f->Power(c,j,&mtrx[j][i]);
 // 	}
 //      delete c;
 //      c=0;
 //    }
 
 //  // print the matrix
 //  //f->MatrixPrint(mtrx,r,r);
  
 //  //////  invert the matrix////////////////////
 //  Element*** invMatrix = 0;
 //  f->MatrixInvert(mtrx,r,&invMatrix);
 //  //cout<<"i am here";
 //  //f->MatrixPrint(invMatrix,r,r);
 //  //cout<<"i am here";

 //  //
 //  Element*** chk = 0;
 //  int rr;
 //  int cc;
 //  //r = 3;
 //  //cout<<"i am here";
 //  f->MatrixMultiply(mtrx,r,r,invMatrix,r,r,&chk,rr,cc);
 //  //f->MatrixPrint(chk,rr,cc);
 //  /////////////////////////////////////////////
 //  Element*** dualBasis=0;
 //  f->MatrixAllocate(&dualBasis,r,1);
 //  //
 //  for(int i = 0; i<r;i++)
 //    {
 //      dualBasis[i][0] = new Element(invMatrix[0][i]);
 //    }
 //  //
 //  unsigned char** bMatrix = 0;
 //  f->ConvertToLowerFieldMatrix(dualBasis,r,1,&bMatrix);
 //  //
 //  unsigned char** bInv = 0;
 //  f->BinaryMatrixInvert(bMatrix,r,&bInv);
 //  //
 //  // cout<<"\n\nbinary matrix\n\n";
 //  // for(int i = 0; i< r; i++)
 //  //   {
 //  //     for(int j = 0;j< r ;j++)
 //  // 	{
 //  // 	  int bv = bInv[i][j];
 //  // 	  cout<<bv<<"    ";
 //  // 	}
 //  //     cout<<endl;
  //   }
  /////////////////////////////////////////////
  // b0 b1 b2
  //int bitToBeMatched = 2;
  //Element * b0 = new Element(f->_zero);
  Element * b2 = 0;//new Element(dualBasis[bitToBeMatched][0]);
  //Element* b1 = new Element(f->_zero);
  //

  for(int allK2Values = 0; allK2Values < r;allK2Values++)
    {
      noOfEqns++;
      //Select c /////////////////////////////////
      Element* k2Value = new Element(f->_dualBasis[allK2Values][0]);
      /////////////////////////////////////////
      //select b2
      long long int val = (long long int)(pow(2,k1Power)+aPower)  % (long long int)(pow(2,fl) - 1);

      long long int kkPower = 0;
      long long int apv = (int)(pow(2,k2Power)+aPower)  % (int)(pow(2,fl) - 1);
      while(val != apv)
	{
	  apv = (long long int)(apv * 2) % (long long int)(pow(2,fl)-1);
	  kkPower = kkPower + 1;//f->_prime;
	}
     kkPower = pow(2,kkPower); 
	//cout<<"\n c power ====>>>  \n"<<cPower;         
      f->Power(k2Value,kkPower,&b2);
      //////////////////////////////////////////
      Element***bCol = 0;//new Element**[3];
      f->MatrixAllocate(&bCol,r,1);
      //
      for(int bi  = 0 ; bi < r;bi++)
	{
	  if(bi == k1Power)
	    {
	      bCol[bi][0] = new Element(b2);
	    }
	  else
	    {
	      if(bi == k2Power)
		{
		  bCol[bi][0]= new Element(k2Value);
		}
	      else
		{
		  bCol[bi][0] = new Element(f->_zero);
		}
	    }
	}
      //
      Element *** tMatrix = 0;
      int tr = 0;
      int tc = 0;
      f->MatrixMultiply(f->_mtrx,r,r,bCol,r,1,&tMatrix,tr,tc);
      unsigned char** tBinaryMatrix = 0;
      f->ConvertToLowerFieldMatrix(tMatrix,tr,tc,&tBinaryMatrix);
      // cout<<"\n\n t binary cofficient matrix\n\n";
      // for(int i = 0; i< r; i++)
      //   {
      //     for(int j = 0;j< r;j++)
      // 	{
      // 	  int bv = tBinaryMatrix[i][j];
      // 	  cout<<bv<<"    ";
      // 	}
      //     cout<<endl;
      //   }
      
      ///////////////////////////////
      // to multiply tI
      unsigned char** aMatrix = 0;
      f->BinaryMatrixMultiply(tBinaryMatrix,r,r,f->_bInv,r,r,&aMatrix);
      cout<<"\n\nbinary cofficient matrix aij\n\n";
      for(int i = 0; i< r; i++)
	{
	  for(int j = 0;j< r;j++)
	    {
	      int bv = aMatrix[i][j];
	      cout<<bv<<"    ";
	    }
	  cout<<endl;
	}
      
      // cout<<"\nbinary coff matrix vj\n\n";
      // for(int i = 0 ; i < r; i++)
      // 	{
      // 	  int bv = vValue[0][i];
      // 	  cout<<bv<<"  ";
      // 	}
      // cout<<endl<<endl;
      // //
      //////////////////////////////////////////////////
      unsigned char s1 = 0;
      unsigned char s2 =0;
      //
      for(int p = 0 ; p < -1;p++)//pow(2,fl)-1;p++)
	{
	  Element* ap = 0;
	  f->Power(f->_a,p,&ap);
	  //
	  Element* tmpInv = 0;
	  if(aPower == -1)
	    {
	      f->Invert(ap,&tmpInv);
	    }
	  else
	    {
	      f->Power(ap,aPower,&tmpInv);
	    }
	  //
	  s1 = 0;
	  for(int i = 0;i<r;i++)
	    {
	      for(int j = 0; j<r;j++)
		{
		  s1 = (~s1 & (aMatrix[i][j] & tmpInv->_val[j] & ap->_val[i]))
		    |
		    (s1 & ~(aMatrix[i][j] & tmpInv->_val[j] & ap->_val[i]));
		  s1 &=1;
		}
	    }
	  //
	  // s2 = 0;
	  // for( int yj = 0; yj < r;yj++)
	  //   {
	  //     s2 =  s2 & ~(vValue[0][yj] & tmpInv->_val[yj]) |
	  // 	(~s2 & (vValue[0][yj] & tmpInv->_val[yj]));
	  //     s2 &= 1;
	  //   }
	  //
	  unsigned char s = (int) s1;//(s1&~s2) | (~s1&s2);
	  s &= 1; 
	  delete ap;
	  ap = 0;
	  delete tmpInv;
	  tmpInv = 0;
	  /////////////////////////////////
	  int v1 = 0;
	  int v2 = (int) s;
	  //cout<<"\n\nv2 = "<<v2<<endl;
	  if(v1 == v2)
	    {
	      //cout<<"\n\n MATCHING.................\n\n";
	    }
	  else{
	    cout<<"\n\nits not working\n\n";
	  }
	  
	}
      /////////////Delete values used
      f->MatrixDelete(&bCol,r,1);
      f->MatrixDelete(&tMatrix,tr,tc);
      //
      delete b2;
      b2 = 0;
      //
      delete k2Value;
      k2Value = 0;
      //
      for(int i = 0;i<r;i++)
	{
	  delete [] aMatrix[i];
	  aMatrix[i] = 0;
	  //
	  delete [] tBinaryMatrix[i];
	  tBinaryMatrix[i] = 0;
	  //
	  //delete [] cBinaryMatrix[i];
	  //cBinaryMatrix[i] = 0;
	}
      
      delete [] aMatrix;
      aMatrix = 0;
      delete [] tBinaryMatrix;
      tBinaryMatrix = 0; 
    }
  //**********************************************************************//
  ///////////////////////////
  // delete the matrix
  // f->MatrixDelete(&mtrx,r,r);
  // f->MatrixDelete(&invMatrix,r,r);
  // f->MatrixDelete(&chk,rr,cc);
  // f->MatrixDelete(&dualBasis,r,1);
  
  //
  // for(int i = 0;i<r;i++)
  //   {
  //     delete [] bMatrix[i];
  //     bMatrix[i] = 0;
  //     //
  //     delete [] bInv[i];
  //     bInv[i] = 0;
  //   }
  // delete [] bMatrix;
  // bMatrix = 0;
  // delete [] bInv;
  // bInv = 0;
  return 1;

}


int SameBigCosetForC(Field* f, int kPower,int aPower)
{
  fflush(stdout);
  int fl = f->_len-1;
  //Field* f = new Field(2,fl+1,(unsigned char*) "1101");// ir = 1+x+x^3
  //Element* a = new Element((unsigned char*)"010",fl);
  //
  //Forming the matrix.....
  int r = fl;
  // Element*** mtrx  =0;// new Element**[r];
 //  f->MatrixAllocate(&mtrx,r,r);
 // //
 //  long int p =1;
 //  for(int i = 0; i< r;i++)
 //    {
 //      Element*c = 0;
 //      f->Power(f->_a,p,&c);
 //      p *= 2;
 //      for(int j = 0 ; j < r;j++)
 // 	{
 // 	  f->Power(c,j,&mtrx[j][i]);
 // 	}
 //      delete c;
 //      c=0;
 //    }
 
 //  // print the matrix
 //  //f->MatrixPrint(mtrx,r,r);
  
 //  //////  invert the matrix////////////////////
 //  Element*** invMatrix = 0;
 //  f->MatrixInvert(mtrx,r,&invMatrix);
 //  //cout<<"i am here";
 //  //f->MatrixPrint(invMatrix,r,r);
 //  //cout<<"i am here";

 //  //
 //  Element*** chk = 0;
 //  int rr;
 //  int cc;
 //  //r = 3;
 //  //cout<<"i am here";
 //  f->MatrixMultiply(mtrx,r,r,invMatrix,r,r,&chk,rr,cc);
 //  //f->MatrixPrint(chk,rr,cc);
 //  /////////////////////////////////////////////
 //  Element*** dualBasis=0;
 //  f->MatrixAllocate(&dualBasis,r,1);
 //  //
 //  for(int i = 0; i<r;i++)
 //    {
 //      dualBasis[i][0] = new Element(invMatrix[0][i]);
 //    }
 //  //
 //  unsigned char** bMatrix = 0;
 //  f->ConvertToLowerFieldMatrix(dualBasis,r,1,&bMatrix);
 //  //
 //  unsigned char** bInv = 0;
 //  f->BinaryMatrixInvert(bMatrix,r,&bInv);
 //  //
 //  // cout<<"\n\nbinary matrix\n\n";
 //  // for(int i = 0; i< r; i++)
 //  //   {
 //  //     for(int j = 0;j< r ;j++)
 //  // 	{
 //  // 	  int bv = bInv[i][j];
 //  // 	  cout<<bv<<"    ";
 //  // 	}
 //  //     cout<<endl;
 //  //   }
  /////////////////////////////////////////////
  // b0 b1 b2
  //int bitToBeMatched = 2;
  //Element * b0 = new Element(f->_zero);
  Element * b2 = 0;//new Element(dualBasis[bitToBeMatched][0]);
  //Element* b1 = new Element(f->_zero);
  //

  for(int allCValues = 0; allCValues < r;allCValues++)
    {
      noOfEqns++;
      //Select c /////////////////////////////////
      Element* cValue = new Element(f->_dualBasis[allCValues][0]);
      /////////////////////////////////////////
      //select b2
      long long int val = (long long int)(pow(2,kPower) + aPower) % (long long int)(pow(2,fl) - 1);
      long long int cPower = 0;
      long long int apv = (aPower == -1)? (long long int)(pow(2,fl) -1 + aPower):aPower;
      //cout<<"\n---ap>"<<aPower<<endl;
      while(val != apv)
	{
	  apv = (long long int)(apv * 2) % (long long int)(pow(2,fl)-1);
	  cPower = cPower +1;// f->_prime;
	}
      cPower = pow(2,cPower);
      //cout<<"\n c power ====>>>  \n"<<cPower;         
      f->Power(cValue,cPower,&b2);
      //delete cValue;
      //cValue = 0;
      // get vj's thereof
      unsigned char** cBinaryMatrix = 0;
      //cout<<"\n\nyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy\n\n";
      //
      f->ConvertToLowerFieldMatrix(cValue,true,&cBinaryMatrix);
      delete cValue ;cValue = 0;
      //cout<<"\n\nyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy\n\n";
      unsigned char** vValue = 0;
      f->BinaryMatrixMultiply(cBinaryMatrix,1,r,f->_bInv,r,r,&vValue);
      //
      //cout<<"\n\n v Matrix \n\n";
      for(int i = 0 ; i < r ; i++)
	{  
//	  int vvv = (int)vValue[0][i];
//	  cout<<vvv<<"     ";
	}
 //     cout<<endl<<endl;
      //
      //////////////////////////////////////////
      Element***bCol = 0;//new Element**[3];
      f->MatrixAllocate(&bCol,r,1);
      //
      for(int bi  = 0 ; bi < r;bi++)
	{
	  if(bi == kPower)
	    {
	      bCol[bi][0] = new Element(b2);
	    }
	  else
	    {
	      bCol[bi][0] = new Element(f->_zero);
	    }
	  //bCol[1][0] = new Element(b1);
	  //bCol[2][0] = new Element(b2);
	}
      //
      Element *** tMatrix = 0;
      int tr = 0;
      int tc = 0;
      f->MatrixMultiply(f->_mtrx,r,r,bCol,r,1,&tMatrix,tr,tc);
      unsigned char** tBinaryMatrix = 0;
      f->ConvertToLowerFieldMatrix(tMatrix,tr,tc,&tBinaryMatrix);
      // cout<<"\n\n t binary cofficient matrix\n\n";
      // for(int i = 0; i< r; i++)
      //   {
      //     for(int j = 0;j< r;j++)
      // 	{
      // 	  int bv = tBinaryMatrix[i][j];
      // 	  cout<<bv<<"    ";
      // 	}
      //     cout<<endl;
      //   }
      
      ///////////////////////////////
      // to multiply tI
      unsigned char** aMatrix = 0;
      f->BinaryMatrixMultiply(tBinaryMatrix,r,r,f->_bInv,r,r,&aMatrix);
      cout<<"\n\nbinary cofficient matrix aij\n\n";
      for(int i = 0; i< r; i++)
	{
	  for(int j = 0;j< r;j++)
	    {
	      int bv = aMatrix[i][j];
	      cout<<bv<<"    ";
	    }
	  cout<<endl;
	}
      
      cout<<"\nbinary coff matrix vj\n\n";
      for(int i = 0 ; i < r; i++)
	{
	  int bv = vValue[0][i];
	  cout<<bv<<"  ";
	}
      cout<<endl<<endl;
      //
      //////////////////////////////////////////////////
      unsigned char s1 = 0;
      unsigned char s2 =0;
      //
      for(int p = 0 ; p < -1;p++)//pow(2,fl)-1;p++)
	{
	  Element* ap = 0;
	  f->Power(f->_a,p,&ap);
	  //
	  Element* tmpInv = 0;
	  if(aPower == -1)
	    {
	      f->Invert(ap,&tmpInv);
	    }
	  else
	    {
	      f->Power(ap,aPower,&tmpInv);
	    }
	  //
	  s1 = 0;
	  for(int i = 0;i<r;i++)
	    {
	      for(int j = 0; j<r;j++)
		{
		  s1 = (~s1 & (aMatrix[i][j] & tmpInv->_val[j] & ap->_val[i]))
		    |
		    (s1 & ~(aMatrix[i][j] & tmpInv->_val[j] & ap->_val[i]));
		  s1 &=1;
		}
	    }
	  //
	  s2 = 0;
	  for( int yj = 0; yj < r;yj++)
	    {
	      s2 =  s2 & ~(vValue[0][yj] & tmpInv->_val[yj]) |
		(~s2 & (vValue[0][yj] & tmpInv->_val[yj]));
	      s2 &= 1;
	    }
	  //
	  unsigned char s = (s1&~s2) | (~s1&s2);
	  s &= 1; 
	  delete ap;
	  ap = 0;
	  delete tmpInv;
	  tmpInv = 0;
	  /////////////////////////////////
	  int v1 = 0;
	  int v2 = (int) s;
	  //cout<<"\n\nv2 = "<<v2<<endl;
	  if(v1 == v2)
	    {
	      //cout<<"\n\n MATCHING.................\n\n";
	    }
	  else{
	    cout<<"\n\nits not working\n\n";
	  }
	  
	}
      /////////////Delete values used
      f->MatrixDelete(&bCol,r,1);
      f->MatrixDelete(&tMatrix,tr,tc);
      //
      delete b2;
      b2 = 0;
      //
      for(int i = 0;i<r;i++)
	{
	  delete [] aMatrix[i];
	  aMatrix[i] = 0;
	  //
	  delete [] tBinaryMatrix[i];
	  tBinaryMatrix[i] = 0;
	  //
	}
      
      delete [] vValue[0];
      vValue[0] = 0;
      delete [] vValue;
      vValue = 0;
      delete [] aMatrix;
      aMatrix = 0;
      delete [] tBinaryMatrix;
      tBinaryMatrix = 0;
      delete [] cBinaryMatrix[0];
      cBinaryMatrix[0] = 0;
      delete [] cBinaryMatrix;
      cBinaryMatrix = 0;
 
    }
  //**********************************************************************//
  ///////////////////////////
  // delete the matrix
  // f->MatrixDelete(&mtrx,r,r);
  // f->MatrixDelete(&invMatrix,r,r);
  // f->MatrixDelete(&chk,rr,cc);
  // f->MatrixDelete(&dualBasis,r,1);
  
  //
  // for(int i = 0;i<r;i++)
  //   {
  //     delete [] bMatrix[i];
  //     bMatrix[i] = 0;
  //     //
  //     delete [] bInv[i];
  //     bInv[i] = 0;
  //   }
  // delete [] bMatrix;
  // bMatrix = 0;
  // delete [] bInv;
  // bInv = 0;
  return 1;
}




int HammingWtOneForC(Field* f,int kPower,int aPower)
{

fflush(stdout);
  int fl = f->_len-1;
  //Field* f = new Field(2,fl+1,(unsigned char*) "1101");// ir = 1+x+x^3
  //Element* a = new Element((unsigned char*)"010",fl);
  //
  //Forming the matrix.....
  int r = fl;
  // Element*** mtrx  =0;// new Element**[r];
 //  f->MatrixAllocate(&mtrx,r,r);
 // //
 //  long int p =1;
 //  for(int i = 0; i< r;i++)
 //    {
 //      Element*c = 0;
 //      f->Power(f->_a,p,&c);
 //      p *= 2;
 //      for(int j = 0 ; j < r;j++)
 // 	{
 // 	  f->Power(c,j,&mtrx[j][i]);
 // 	}
 //      delete c;
 //      c=0;
 //    }
 
 //  // print the matrix
 //  //f->MatrixPrint(mtrx,r,r);
  
 //  //////  invert the matrix////////////////////
 //  Element*** invMatrix = 0;
 //  f->MatrixInvert(mtrx,r,&invMatrix);
 //  //f->MatrixPrint(invMatrix,r,r);
 //  //
 //  Element*** chk = 0;
 //  int rr;
 //  int cc;
 //  f->MatrixMultiply(mtrx,r,r,invMatrix,r,r,&chk,rr,cc);
 //  //f->MatrixPrint(chk,rr,cc);
 //  /////////////////////////////////////////////
 //  Element*** dualBasis=0;
 //  f->MatrixAllocate(&dualBasis,r,1);
 //  //
 //  for(int i = 0; i<r;i++)
 //    {
 //      dualBasis[i][0] = new Element(invMatrix[0][i]);
 //    }
 //  //
 //  unsigned char** bMatrix = 0;
 //  f->ConvertToLowerFieldMatrix(dualBasis,r,1,&bMatrix);
 //  //
 //  unsigned char** bInv = 0;
 //  f->BinaryMatrixInvert(bMatrix,r,&bInv);
 //  //
 //  // cout<<"\n\nbinary matrix\n\n";
 //  // for(int i = 0; i< r; i++)
 //  //   {
 //  //     for(int j = 0;j< r ;j++)
 //  // 	{
 //  // 	  int bv = bInv[i][j];
 //  // 	  cout<<bv<<"    ";
 //  // 	}
 //  //     cout<<endl;
 //  //   }
 //  /////////////////////////////////////////////
  // b0 b1 b2
  

  //Select c /////////////////////////////////
  for(int cIndex = 0; cIndex < r;cIndex++)
    {
      noOfEqns++;
      //cout<<"\n cIndex ============="<<cIndex<<endl;
      Element* cValue = new Element(f->_dualBasis[cIndex][0]);
      // get vj's thereof
      unsigned char** cBinaryMatrix = 0;
      f->ConvertToLowerFieldMatrix(cValue,true,&cBinaryMatrix);
  
      unsigned char** vValue = 0;
      f->BinaryMatrixMultiply(cBinaryMatrix,1,r,f->_bInv,r,r,&vValue);
      //
      cout<<"\n\n v cofficient matrix Matrix vj \n\n";
      for(int i = 0 ; i < r ; i++)
	{  
	  int vvv = (int)vValue[0][i];
	  cout<<vvv<<"     ";
	}
      cout<<endl<<endl;
      //
//TO DO: write the rhs
      unsigned char* rhsBit = new unsigned char[f->_len-1];
      cout<<"\n The rightHand side xbis  matrix :\n";
      for(int rhs = 0;rhs < fl;rhs++)
	{
	  //
	  Element* ap = 0;
	  f->Power(f->_a,rhs,&ap);
	  ///////////////////////////////////////////////////////
	  Element* traceOf=0;
	  Element*traceValue = 0;
	  Element* xPower = 0;
	  int pr = (int)(aPower) % (int)(pow(2,fl) -1);
	  f->Power(ap,pr,&xPower);
	  f->Multiply(cValue,xPower,&traceOf);
	  f->Trace(traceOf,&traceValue);
	  unsigned char xBit = traceValue->_val[0];
	  delete traceOf,traceValue,xPower;
	  traceOf = 0;traceValue = 0;xPower = 0;
	  rhsBit[rhs]=xBit;
	  int vv = xBit;
	  cout<<vv<<"   ";
	}
      cout<<"\n\n";
      delete []rhsBit;rhsBit = 0;
      
      //////////////////////////////////CHECK/////////////////////////////////////
      unsigned char bit = 0;
      for(int i = 0; i< -1;i++)//pow(2,fl)-1;i++)
	{
	  Element* ap = 0;
	  Element* apy = 0;
	  f->Power(f->_a,i,&ap);
	  f->Power(ap,aPower,&apy);
	  Element* traceOf =0;
	  f->Multiply(cValue,apy,&traceOf);
	  Element* traceValue = 0;
	  f->Trace(traceOf,&traceValue);
	  bit = traceValue->_val[0];
	  //
	  unsigned char s2 = 0;
	  for( int yj = 0; yj < r;yj++)
	    {
	      s2 =  s2 & ~(vValue[0][yj] & apy->_val[yj]) |
		(~s2 & (vValue[0][yj] & apy->_val[yj]));
	      s2 &= 1;
	    }
	  int v1 = (int)bit;
	  int v2 = (int)s2;
	  if(v1 != v2)
	    {
	      cout<<"\n\n its not matching \n\n";
	    }    
	  //
	  delete ap,apy,traceOf,traceValue;
	  ap=0;apy=0;traceOf = 0;traceValue = 0;
	}

      ////////////delete section//////////////////////////////
      delete cValue ;cValue = 0;
      
      //for(int i = 0;i<r;i++)
      //{
	  delete[] cBinaryMatrix[0];
	  cBinaryMatrix[0] = 0;
	  //}
      delete [] cBinaryMatrix;
      cBinaryMatrix = 0;
      //
      delete [] vValue[0];
      vValue[0] = 0;
      delete [] vValue;
      vValue = 0;
      //cout<<"\n*******************************************\n";
    }
  //////////////////////////////////////////////////
  // delete the matrix
  // f->MatrixDelete(&mtrx,r,r);
  // f->MatrixDelete(&invMatrix,r,r);
  // f->MatrixDelete(&chk,rr,cc);
  // f->MatrixDelete(&dualBasis,r,1);

  //
  // for(int i = 0;i<r;i++)
  //   {
  //     // delete [] bMatrix[i];
  //     //bMatrix[i] = 0;

  //     //delete [] cBinaryMatrix[i];
  //     //cBinaryMatrix[i] = 0;
     
  //     //
  //     //delete [] bInv[i];
  //     //bInv[i] = 0;
  //     //
  //   }

  // delete [] bMatrix;
  // bMatrix = 0;
  // delete [] bInv;
  //bInv = 0;
  return 1;
}

///////////////////////////////////////////
int HammingWtOneForBk(Field* f , int kPower,int aPower)
{
  fflush(stdout);
  int fl = f->_len-1;
  //Field* f = new Field(2,fl+1,(unsigned char*) "11001");// ir = 1+x+x^3
  //Element* a = new Element((unsigned char*)"0100",fl);
  //
  //Forming the matrix.....
  int r = fl;
  // Element*** mtrx  =0;// new Element**[r];
  // f->MatrixAllocate(&mtrx,r,r);
  // //
  // long int p =1;
  // for(int i = 0; i< r;i++)
  //   {
  //     Element*c = 0;
  //     f->Power(f->_a,p,&c);
  //     p *= 2;
  //     for(int j = 0 ; j < r;j++)
  // 	{
  // 	  f->Power(c,j,&mtrx[j][i]);
  // 	}
  //     delete c;
  //     c=0;
  //   }
 
  // // print the matrix
  // //f->MatrixPrint(mtrx,r,r);
  
  // //////  invert the matrix////////////////////
  // Element*** invMatrix = 0;
  // f->MatrixInvert(mtrx,r,&invMatrix);
  // //cout<<"i am here";
  // //f->MatrixPrint(invMatrix,r,r);
  // //cout<<"i am here";

  // //
  // Element*** chk = 0;
  // int rr;
  // int cc;
  // //cout<<"i am here";
  // f->MatrixMultiply(mtrx,r,r,invMatrix,r,r,&chk,rr,cc);
  // //f->MatrixPrint(chk,rr,cc);
  // /////////////////////////////////////////////
  // Element*** dualBasis=0;
  // f->MatrixAllocate(&dualBasis,r,1);
  // //
  // for(int i = 0; i<r;i++)
  //   {
  //     dualBasis[i][0] = new Element(invMatrix[0][i]);
  //   }
  // //
  // unsigned char** bMatrix = 0;
  // f->ConvertToLowerFieldMatrix(dualBasis,r,1,&bMatrix);
  // //
  // unsigned char** bInv = 0;
  // f->BinaryMatrixInvert(bMatrix,r,&bInv);

  // // unsigned char** iM = 0;
  // // f->BinaryMatrixMultiply(bMatrix,r,r,bInv,r,r,&iM);
  // // //
  // // cout<<"\n\nbinary matrix\n\n";
  // // for(int i = 0; i< r; i++)
  // //   {
  // //     for(int j = 0;j< r ;j++)
  // // 	{
  // // 	  int bv = (int)iM[i][j];
  // // 	  cout<<bv<<"    ";
  // // 	}
  // //     cout<<endl;
  // //   }
  // // //
  // // for(int i = 0;i<r;i++)
  // //   {
  // //     delete [] iM[i];
  // //   }
  // // delete [] iM;
  // // iM = 0;
  
  /////////////////////////////////////////////
  // b0 b1 b2
  
  for(int bitToBeMatched = 0;bitToBeMatched < r;bitToBeMatched++)
    {
      noOfEqns++;
      Element***bCol = 0;
      f->MatrixAllocate(&bCol,r,1);
      //
      for(int bi = 0 ; bi < r;bi++)
	{
	  if(bi == kPower)
	    {
	      bCol[bi][0] = new Element(f->_dualBasis[bitToBeMatched][0]);
	    }
	  //
	  else
	    {
	      bCol[bi][0] = new Element(f->_zero);
	    }
      
	}
  //
      Element *** tMatrix = 0;
      int tr = 0;
      int tc = 0;
      f->MatrixMultiply(f->_mtrx,r,r,bCol,r,1,&tMatrix,tr,tc);
      unsigned char** tBinaryMatrix = 0;
      f->ConvertToLowerFieldMatrix(tMatrix,tr,tc,&tBinaryMatrix);
      //
      // cout<<"\n\n t binary cofficient matrix\n\n";
      // for(int i = 0; i< r; i++)
      //   {
      //     for(int j = 0;j< r;j++)
      // 	{
      // 	  int bv = tBinaryMatrix[i][j];
      // 	  cout<<bv<<"    ";
      // 	}
      //     cout<<endl;
      //   }
      
      ///////////////////////////////
      // to multiply tI
      unsigned char** aMatrix = 0;
      f->BinaryMatrixMultiply(tBinaryMatrix,r,r,f->_bInv,r,r,&aMatrix);
      cout<<"\n\nbinary cofficient matrix aij\n\n";
      for(int i = 0; i< r; i++)
	{
	  for(int j = 0;j< r;j++)
	    {
	      int bv = aMatrix[i][j];
	      cout<<bv<<"    ";
	    }
	  cout<<endl;
	}
      //
      //TO DO: write the rhs
      unsigned char* rhsBit = new unsigned char[f->_len-1];
      cout<<"\n The rightHand side matrix :\n";
      for(int rhs = 0;rhs < fl;rhs++)
	{
	  //
	  Element* ap = 0;
	  f->Power(f->_a,rhs,&ap);
	  ///////////////////////////////////////////////////////
	  Element* traceOf=0;
	  Element*traceValue = 0;
	  Element* xPower = 0;
	  int pr = (int)(pow(2,kPower)+aPower) % (int)(pow(2,fl) -1);
	  f->Power(ap,pr,&xPower);
	  f->Multiply(bCol[kPower][0],xPower,&traceOf);
	  f->Trace(traceOf,&traceValue);
	  unsigned char xBit = traceValue->_val[0];
	  delete traceOf,traceValue,xPower;
	  traceOf = 0;traceValue = 0;xPower = 0;
	  rhsBit[rhs]=xBit;
	  int vv = xBit;
	  cout<<vv<<"   ";
	}
      cout<<"\n\n";
      delete []rhsBit;rhsBit = 0;
      //////////////CHECK////////////////////////////////////
      unsigned char s = 0;
      unsigned char xBit = 0;
      for(int p = 0 ; p < -1;p++)//pow(2,fl)-1;p++)
	{
	  Element* ap = 0;
	  f->Power(f->_a,p,&ap);
	  Element* apy=0;
	  if(aPower == -1)
	    {
	      f->Invert(ap,&apy);
	    }
	  else
	    {
	      f->Power(ap,aPower,&apy);
	    }
	  ///////////////////////////////////////////////////////
	  
	  Element* traceOf=0;
	  Element*traceValue = 0;
	  Element* xPower = 0;
	  int pr = (int)(pow(2,kPower)+aPower) % (int)(pow(2,fl) -1);
	  f->Power(ap,pr,&xPower);
	  f->Multiply(bCol[kPower][0],xPower,&traceOf);
	  f->Trace(traceOf,&traceValue);
	  xBit = traceValue->_val[0];
	  delete traceOf,traceValue,xPower;
	  traceOf = 0;traceValue = 0;xPower = 0;
	  ///////////////////////////////////////////////////////
	  s = 0;
	  for(int i = 0; i< r;i++)
	    {
	      for(int j = 0; j<r;j++)
		{
		  s =  (s & ~( aMatrix[i][j] & ap->_val[i] & apy->_val[j])) | (~s &(aMatrix[i][j] & ap->_val[i] & apy->_val[j]));
		}
	    }
	  
	  int v1 = (int)xBit;
	  int v2 = (int)s;
	  if(v1 == v2)
	    {
	      //cout<<"\n\n MATCHING.................\n\n";
	    }
	  else{
	    ap->Show();cout<<endl;apy->Show();
	    cout<<"\n\n not matching for p ="<<p<<"\n\n";
	  }
	  delete ap;
	  ap = 0;
	  delete apy;
	  apy=0;
	}
      ////////////////////////
      // delete parameters used in loop
      f->MatrixDelete(&bCol,r,1);
      for(int i = 0;i<r;i++)
    {
      delete [] aMatrix[i];
      aMatrix[i] = 0;
      //
      delete [] tBinaryMatrix[i];
      tBinaryMatrix[i] = 0;
    }
      delete [] aMatrix;
      aMatrix = 0;
      delete [] tBinaryMatrix;
      tBinaryMatrix = 0;
      f->MatrixDelete(&tMatrix,tr,tc);
    }
  //********************************************************************//
  //////////////////////////////////////////////////
  // delete the matrix
  // f->MatrixDelete(&mtrx,r,r);
  // f->MatrixDelete(&invMatrix,r,r);
  // f->MatrixDelete(&chk,rr,cc);
  // f->MatrixDelete(&dualBasis,r,1);
  
 
  //
  // for(int i = 0;i<r;i++)
  //   {
  //     //
  //     delete [] bMatrix[i];
  //     bMatrix[i] = 0;
  //     //
  //     delete [] bInv[i];
  //     bInv[i] = 0;
         
  //   }

  // delete [] bMatrix;
  // bMatrix = 0;
  // delete [] bInv;
  //bInv = 0;
  return 1;
}


int inv_main()
{
  fflush(stdout);
  Field *f = new Field(2,9,(unsigned char*) "101110001");// ir = 1+x^2+x^5
  //
  Element * a = new Element((unsigned char*)"00111010",(size_t)8);
  Element* b = 0;
  f->Invert(a,&b);
  
cout<<"\ninverse======================>>>>\n";

  b->Show();
  cout<<endl;
  Element* c = 0;
  f->Multiply(a,b,&c);
  //
  cout<<"\n---chkkkk--\n" ;c->Show();delete c; c = 0;
  delete a;
  a = 0;
  delete b;
  b=0;
  return 0;
}



int GenerateField()
{
  fflush(stdout);
  Field *f = new Field(2,9,(unsigned char*) "101110001");// ir = 1+x^2+x^5
  //
  Element * a = new Element((unsigned char*)"01000000",(size_t)8);
  Element * b = new Element((unsigned char*)"10000000",(size_t)8);
  Element*c = 0;
  Element*d = 0;
  Element*e = 0;
  Element* g =0;
  Element* h=0;
  Element** aryEls =  new Element* [255];
  //
  for(int i = 0; i< 255;i++)
    {
      aryEls[i] = 0;
    }
  //
  cout<<" power           alpha               f(x)              F(x+1)               f(x)+f(x+1)"<<endl;
  cout<<"-----------------------------------------------------------------------------------------------"<<endl;
  for(int i = 1;i<=255;i++)
    {
      f->Power(a,i,&c);
      
      f->Power(c,39,&d);
      cout<<" "<<i << "\t| "; 
      c->Show();     //////////////// alpha ^ i
      cout<<"\t|\t ";
      d->Show();     ///////////////  alpha ^ i ^ a
      //
      f->Add(c,b,&e);
      f->Power(e,39,&g);
      //
      cout<< "\t|\t ";
      g->Show();              
      //
      cout<< "\t|\t ";
      f->Add(d,g,&h);
      h->Show();
      aryEls[i-1] = new Element(h);
      cout<<endl;
      delete d;
      d = 0;
      delete c;
      c= 0;
      delete e;
      e = 0;
      delete g;
      g = 0;
      delete h;
      h=0;
    }

  //
  for(int i = 0; i <  255;i++)
    {
      for(int j = 254;j>i;j--)
	{
	  if(aryEls[j]== 0 || aryEls[j-1] == 0)
	    {
	      continue;
	    }
	  if(H(aryEls[j]) < H(aryEls[j-1]))
	    {
	      Element* temp = new Element(aryEls[j]);
	      delete aryEls[j];
	      aryEls[j] = 0;
              aryEls[j] = new Element(aryEls[j-1]);
              delete  aryEls[j-1];
	      aryEls[j-1] = 0;
	      aryEls[j-1] = new Element(temp);
	      delete temp;
	      temp = 0;
	    }
	}
    }

  /////////////////////////////////////////////////
  cout<<"\n\n sorted values \n\n";
  for(int i = 0 ; i < 255;i++)
    {
      if(aryEls[i] == 0)
	{
	  continue;
	}
      aryEls[i]->Show();
      cout<<endl;
    }

  ////////////////////////////////////////////////
  delete f;
  delete a;
  delete b;
  delete c;
  delete e;
  delete g;
  for(int i = 0;i< 255;i++)
    {
      if(aryEls[i] != 0)
	{
	  delete aryEls[i];
          aryEls[i] = 0;
	}
    }
  //
  delete [] aryEls;
  aryEls = 0;
}

///////////////////////////////////////////////////////////////////////

class CosetManager
{
public:
  bool _bValue;
  double _kValue;
  double  _cLeader;
  double  _cSize;
  CosetManager(bool b,int v ,int  cl , int cs)
  {
    _bValue = b;
    _kValue = v;
    _cLeader = cl;
    _cSize = cs;

  }
  CosetManager(CosetManager* m)
  {
    _bValue = m->_bValue;
    _kValue = m->_kValue;
    _cLeader = m->_cLeader;
    _cSize = m->_cSize;
  }
  //
  ~CosetManager()
  {
  }
  //
  bool Print()
  {
    long long int d = pow(2,_kValue)+9;
    if(_bValue == false)
      {
	d = _kValue;
      }
    cout<<_cLeader<<"\t"<<_kValue<<"\t"<<_cSize<<"\t"<<d<<endl;
    return true;  
  }
};


int ComputeBiAffineEquations()
{
  fflush(stdout);
  //Field *f = new Field(2,9,(unsigned char*) "101110001");// 
  //Field *f = new Field(2,11,(unsigned char*) "10010000001");// 
  //Field *f = new Field(2,5,(unsigned char*) "11001");// ir = 1+x^2+x^5
  //Field *f = new Field(2,4,(unsigned char*) "1101");// ir = 1+x^2+x^5
  //Field* f = new Field(2,6,(unsigned char*)"101001");
  //Field* f = new Field(2,21,(unsigned char*)"100100000000000000001"); 
 //Field* f = new Field(2,31,(unsigned char*)"1000000001100000000000000001001"); 
  // Field* f = new Field(2,10,(unsigned char*)"1000100001"); 
   Field* f = new Field(2,7,(unsigned char*)"1100001"); 
  int order = pow(2,f->_len-1)-1;
  CosetManager** csm = new CosetManager*[f->_len];
  int a = 5;
  int aPower = 0;
  if(a > 0)
    aPower =  a % order ;
  else
   aPower = pow(2,f->_len-1)-1+a;

  //cout<<"\naPower = "<<aPower<<endl;

  for(int k = 0;k< f->_len-1;k++)
    {
      int xPower = (int)(pow(2,k)+a) % order;
      int cs = f->GetCosetSize(xPower);
      int cl = f->GetCosetLeader(xPower);
      csm[k] = new CosetManager(true,k,cl,cs);
    }  
  csm[f->_len-1] = new CosetManager(false,aPower,f->GetCosetLeader(aPower),f->GetCosetSize(aPower));


  //
  for(int i = 0;i<f->_len;i++)
    {
      for( int j = f->_len-1;j>i;j--)
	{
	  if(csm[j-1]->_cLeader > csm[j]->_cLeader)
	    {
	      CosetManager* tmp = new CosetManager(csm[j-1]);
	      delete csm[j-1];
	      csm[j-1]=0;
	      csm[j-1]= new CosetManager(csm[j]);
	      delete csm[j];csm[j]=0;
	      csm[j] = new CosetManager(tmp);
	      delete tmp;tmp=0;
	    } 
	}
    }
  //********************Main algorithm***************************************//
  int xPower;

  for(int k = 0;k<f->_len;k++)
    {
      if(csm[k]->_bValue)
	{
	  xPower = (int)(pow(2,k)+a) % order;
	}
      else
	{
	  //xPower = 0;
	}
      
      // IF HAMMING WEIGHT IS 1
       if(f->GetHammingWt(csm[k]->_cLeader) == 1)
       	{
	 
	  if(csm[k]->_bValue)
	    {
               cout<<"\nHamming wt 1 for bk\n";
	      int i = HammingWtOneForBk(f,csm[k]->_kValue,aPower);
          
	      //cout<<"\n\nk is ------------> "<<k<<endl;
	    }
	  else
	    {
	      cout<<"\nHamming wt 1 for C\n";
	      int i = HammingWtOneForC(f,k,csm[k]->_kValue);// =aPower
	    }
      	}
       //////HAMMING WT 1 ENDS
       //ELSE
       else
	 {
	   ////SMALLER COSET CASE
	   if(csm[k]->_cSize < f->_len-1)
	     {
	       if(csm[k]->_bValue)
		 {
		   cout<<"\nSmaller coset for bk\n";
		   SmallerCosetForBk(f,csm[k]->_kValue,aPower);
		   cout<<"\n now k is :: "<<k<<endl;
		 }
	       else
		 {
		   cout<<"\nSmaller coset for C\n";
		   SmallerCosetForC(f,k,csm[k]->_kValue);//=aPower
		   //cout<<"\n aPower is "<< aPower<<endl;
		 }
	     }
	   ///////////////SMALLER COSET CASE ENDS
	   
	   //SAME COSET//////////////////////
	   if(k== f->_len-1)
	     {
	       continue;
	     }
	   if(csm[k]->_cLeader == csm[k+1]->_cLeader)
	     {
	       //BIG SAME COSET
	       if(csm[k]->_cSize == f->_len-1)
		 {
		   // BOTH ARE BK
		   if( (csm[k]->_bValue == true) && (csm[k+1]->_bValue == true))
		     {
		       //noOfEqns = noOfEqns + f->_len-1;
                       cout<<"\nSame BIG coset for bk\n";
		       int i = SameBigCosetForBk(f,csm[k]->_kValue,csm[k+1]->_kValue,aPower);
		     }
		   else // ONE BK ONE C
		     {
		       if((csm[k]->_bValue == true) && csm[k+1]->_bValue == false)
			 {
			   cout<<"\nSame BIG coset for C 1\n"<<csm[k+1]->_kValue<<endl;
			   int i = SameBigCosetForC(f,csm[k]->_kValue,csm[k+1]->_kValue);
			 }
		       else
			 {
			   cout<<"\nSame BIG coset for C 2\n";
			   int i = SameBigCosetForC(f,csm[k+1]->_kValue,csm[k]->_kValue); 
			 }
		     }
		 }
	       //////
	       else////Small
		 {
		   cout<<"\n\n\n\n Kishan kishan da's special case\n\n\n\n";
		 }
	     }

	 }
    }
  
  //*******************************************************************//
  


  ////////DELETE SECTION//////////////////////////////////////////////////
  cout<<"\nl\tk\ts\n";
  for(int i = 0 ; i <f->_len;i++)
    {
      csm[i]->Print();
      delete csm[i];csm[i]=0;
    }
  //
  delete [] csm;
  csm = 0;
  //
  delete f;
  f = 0;
  cout<<"\nNUMBER OF EQUATIONS : "<<noOfEqns<<endl<<endl;
  return 1;
  
}


void Detect(Field* f)
{
  fflush(stdout);
  //Field *f = new Field(2,11,(unsigned char*)"10101010111");
  //Field *f = new Field(2,9,(unsigned char*) "101110001");// 
  // Field *f = new Field(2,5,(unsigned char*) "11001");// ir = 1+x^2+x^5
  //Field *f = new Field(2,4,(unsigned char*) "1101");// ir = 1+x^2+x^5
  //
  int order = (int)(pow(2,f->_len-1)-1);
  for(int a = 1; a<pow(2,f->_len-1)-1;a++)
    {
      CosetManager** csm = new CosetManager*[f->_len];
      int aPower = 0;
      if(a > 0)
	aPower =  a % order ;
      else
	aPower = pow(2,f->_len-1)-1+a;
      
      //cout<<"\naPower = "<<aPower<<endl;
      
      for(int k = 0;k< f->_len-1;k++)
	{
	  int xPower = (int)(pow(2,k)+a) % order;
	  int cs = f->GetCosetSize(xPower);
	  int cl = f->GetCosetLeader(xPower);
	  csm[k] = new CosetManager(true,k,cl,cs);
	}  
      csm[f->_len-1] = new CosetManager(false,aPower,f->GetCosetLeader(aPower),f->GetCosetSize(aPower));
        
      //
      for(int i = 0;i<f->_len;i++)
	{
	  for( int j = f->_len-1;j>i;j--)
	    {
	      if(csm[j-1]->_cLeader > csm[j]->_cLeader)
		{
		  CosetManager* tmp = new CosetManager(csm[j-1]);
		  delete csm[j-1];
		  csm[j-1]=0;
		  csm[j-1]= new CosetManager(csm[j]);
		  delete csm[j];csm[j]=0;
		  csm[j] = new CosetManager(tmp);
		  delete tmp;tmp=0;
		} 
	    }
	}
      
      /////////////////////////////////////////
      bool flg = false;    
      for(int i = 1;i<f->_len;i++)
	{
	  if(csm[i]->_cLeader == csm[i-1]->_cLeader)
	    {
	      if(csm[i]->_cSize < f->_len-1)
		{

		  flg = true;
		  cout<<"\n n = "<<f->_len-1<<"\t a = "<<a<<"\tk = "<<csm[i]->_kValue<<endl;
                  scanf("%d",&a);
		}
	    }
	}
      if(!flg)
	{
	  //cout<<"n= "<<f->_len-1<<"a ="<<a <<"\n\nDOES NOT EXIST\n";
	}
      ////////DELETE SECTION//////////////////////////////////////////////////
      //cout<<"\nl\tk\ts\n";
      for(int i = 0 ; i <f->_len;i++)
	{
	  //csm[i]->Print();
	  delete csm[i];csm[i]=0;
	}
      //
      delete [] csm;
      csm = 0;
      //
    }
  //delete f;
  //f = 0;
  
}

int Count_main()
{
  fflush(stdout);
  noOfEqns = 0;
  //int i = pow(2,16);
  for(int i = 3;i<22;i++)
    {
  //Field *f = new Field(2,9,(unsigned char*) "101110001");// 
   Field *f = new Field(2,i,(unsigned char*) "1101");// ir = 1+x^2+x
  // int i = SmallerCosetForBk(f,0,-1);  
  //int i = SameBigCosetForC(f,2,-1);

  //int i = ComputeBiAffineEquations();
    Detect(f);
  delete f;f = 0;
  cout<<"\n n = "<<i-1<<"\n-----------------------------\n";  //
    }
  return 1;
}


int eqn_main()
{
//Field *f = new Field(2,34,(unsigned char*) "1000000000000100000000000000000001");
//  Field* f = new Field(2,21,(unsigned char*)"100100000000000000001"); 
//Field *f = new Field(2,5,(unsigned char*) "11001");
int i = ComputeBiAffineEquations();
  //int count = Count_main();
// delete f;f=0;  
return 1;
}

  
bool IsAllNonZeroMatrix(Field*f, Element***a,int ar,int ac)
{
    bool bRetVal = true;
    for(int oi = 0 ; oi< ar;oi++)
    {
        for(int oj = 0 ; oj < ac;oj++)
        {
            if(f->Equals(a[oi][oj],f->_zero))
            {
                return false;
            }
        }
    }
    //
    return bRetVal;

}


bool GetMinorMatrix(Field* f,int r, int c, Element*** a,int ar,int ac, Element**** aM)
{

    bool bRetVal = false;
    f->MatrixAllocate(aM,ar-1,ac-1);
    //
    for(int i = 0; i< ar;i++)
    {
        if(i == r)
        {
            continue;
        }
        for(int j = 0 ; j < ac;j++)
        {
            if(j == c)
            {
                continue;
            }
            /////////////////////////////////////////
            int aMi = (i < r)?i:i-1;
            int aMj = (j < c)?j:j-1;
            //cout<<"("<<aMi<<","<<aMj<<")      ";
            (*aM)[aMi][aMj] = new Element(a[i][j]);
            /////////////////////////////////////////
        }
        //cout<<endl;
    }
    //cout<<"\n---------------------------------------------------------------------------\n";
    ///////////////////////////////////////////////////////
    bRetVal = true;
    return bRetVal;
}

bool CheckMDS(Field*f, Element***a,int ar,int ac)
{
    bool bRetVal = true;
    Element*** aInv = 0;
    Element*** aMinor =0;
    Element*** aMinorInv =0;
    //f->MatrixAllocate(&aMinor,ar-1,ac-1);
    ////////////////////////////////////////////////////////////////////////////
    if(!IsAllNonZeroMatrix(f,a,ar,ac))
    {
        return false;
    }
    ////////////////////////////////////////////////////////////////////////////
    //
    f->MatrixInvert(a,ar,&aInv);
    if(!IsAllNonZeroMatrix(f,aInv,ar,ac))
    {
        f->MatrixDelete(&aInv,ar,ac);
        aInv = 0;
        return false;
    }
    ////////////////////////////////////////////////////////////////////////////
    //
    for(int i = 0; i < ar; i++)
    {
        for(int j = 0; j< ac;j++)
        {
            if(aMinor != 0)
            {
                f->MatrixDelete(&aMinor,ar-1,ac-1);
                aMinor=0;
            }

            GetMinorMatrix(f,i,j,a,ar,ac,&aMinor);
            
            if(aMinorInv != 0)
            {
                f->MatrixDelete(&aMinorInv,ar-1,ac-1);
                aMinorInv=0;
            }
            f->MatrixInvert(aMinor,ar-1,&aMinorInv);
            if(!IsAllNonZeroMatrix(f,aMinorInv,ar-1,ac-1))
            {
                bRetVal = false;
                break;
            }
           
        }
        //
        if(bRetVal == false)
        {
            break;
        }
    }
    //
    ////////////////////////////////////////////////////////////////////////////
    // DELETE SECTION
    f->MatrixDelete(&aMinor,ar-1,ac-1);
    f->MatrixDelete(&aInv,ar,ac);
    aMinor=0;
    a=0;
    if(aMinorInv != 0)
     {
        f->MatrixDelete(&aMinorInv,ar-1,ac-1);
        aMinorInv=0;
     }
     //////////////////////////////////////////////////////////////////////////
    return bRetVal;

}

int main()
{
  //Field* f = new Field(2,21,(unsigned char*)"100100000000000000001"); 
  Field* f = new Field(2,9,(unsigned char*)"110110001"); 
  Element*** a = 0;
  Element***b = 0;
  int br,bc;
  f->MatrixAllocate(&a,4,4);
  //////////////////////////////////////////////////
  for(int i = 0; i<3;i++)
  {
      for(int j = 0; j < 4;j++)
      {
          if(j == 0)
          {
              a[i][j] = new Element(f->_zero);
          }
          else{
              
              
                  if(i == j-1)
                  {
                      a[i][j] = new Element(f->_one);
                  }
                  else
                  {
                      a[i][j] = new Element(f->_zero);
                  }
              
          }
      
      }
  }
  //////////////////////////////////////////////////////////////
  a[3][0] = new Element(f->_one);
  a[3][1] = new Element(f->_a);
  a[3][2] = new Element(f->_one);
  f->Power(f->_a,11,&a[3][3]);
  //////////////////////////////////////////////////////////////
  //f->MatrixPrint(a,4,4);
  //////////////////////////////////////////////////
  f->MatrixPower(a,4,4,4,&b,br,bc);
  f->MatrixPrint(b,4,4);
  //
  
  bool bdn =  CheckMDS(f,b,4,4);
  if(bdn)
  {
      cout<<"\nthis is MDS matrix\n\n";
  }
  else
  {
      cout<<"\n NOT MDS\n\n";
  }
  ///////////////////////////////////////////
  f->MatrixDelete(&a,4,4);
  f->MatrixDelete(&b,4,4);
  delete f;
  f=0;

}
