///////////////////////////////////////////////////////////////////////////////////////
//!                                                          //////////////////////////
//!Author : Indranil Ghosh Ray...............................//////////////////////////
//!Purpose: Fun..............................................//////////////////////////
//!Subject: Finite field.....................................//////////////////////////
//!Language:c++..............................................//////////////////////////
//!                                                          //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<string>
int noOfEqns;
int noOfHwtEqns;
int noOfSmallerCosetEqns;
int noOfSameCosetEqns;
using namespace std;

class Field;

class Element
{
  friend class Field;
public:
  // c/dTor
  Element(){}
  Element(unsigned char* val,size_t l);
  Element(Element*);
  ~Element();
  //////////////
  bool Show();
  //
  friend int H(Element*);
  //private:
  size_t _len;
  unsigned char* _val; 
};

int H(Element* e)
{
  int t = 0;
  int m = 1;
  for(int i = e->_len-1;i>=0;i--)
    {
      int vl = (int) e->_val[i];
      if(vl>=48)
	{
	  vl = vl - 48;
	}
      t += (m*vl);
      m *= 2;
    }
  return t;
}

Element::Element(Element* e)
{
  //
  _len = e->_len;
  _val = new unsigned char[_len];
  memcpy(_val,e->_val,_len);

}

Element::Element(unsigned char* v,size_t l)
{
  _len = l;
  _val = new unsigned char[l];
  memcpy(_val,v,l);
   for(int i = 0 ; i < l ; i++)
     {
       _val[i] = _val[i] &1;
     }
}

Element::~Element()
{
  if(_val !=0)
    {
      delete [] _val;
      _val = 0;
    }
}

bool Element::Show()
{
  bool bRetVal = false;
  if(_val == 0)
    {
      cout<<"\n\noops! this is void\n";
      return bRetVal;
    }
  for(int i = 0; i< _len ;i++)
    {
      //_val[i] &= 1;
      int v = _val[i]&1;
      v = (v>=48)?v-48:v; 
      cout<<v;
    }
}
/// Class "Field" declaration ///////////////////////////////////////
class Field
{
  //friend class Element;
public:
  // ! C/DTor//////////
  Field(){}
  Field(int prime,size_t len ,const unsigned char* ir);
  ~Field();
  bool LoadBasis();
  bool FreeBasis();
  //
  bool Add(Element* a, Element* b,Element** c);
  bool Multiply(Element* a, Element* b,Element** c,bool d);
  bool Power(Element *a,long long int p,Element** b);
  bool PolynomialAdd(Element* a,Element* b, Element** out);
  bool PolynomialDivide(const Element* el ,const Element*v,Element**u,Element** r);
  bool Invert(Element* a,Element** out);
  bool Equals(Element* a,Element* b);
  bool MatrixInvert(Element*** a,int r,Element**** out);
  bool MatrixDelete(Element****b,int r,int c); 
  bool MatrixPrint(Element*** mtrx,int r,int c);
  bool MatrixAllocate(Element**** b,int r,int c);
  bool MatrixMultiply
  (
   Element***a, int ar,int ac,
   Element***b,  int br,int bc,
   Element****c, int&cr,int&cc
   );
  bool MatrixPower(Element***a,int ar,int ac,int pow,Element****b, int& br, int& bc);

  bool Trace(Element*a,Element**b,int m);
  bool ConvertToLowerFieldMatrix(Element***e,int r,int c,unsigned char*** out);
  bool ConvertToLowerFieldMatrix(Element*el,bool rMatrix ,unsigned char*** out);
  bool BinaryMatrixInvert(unsigned char**a,int r,unsigned char*** out);
  bool BinaryMatrixMultiply
(unsigned char**a,int ra,int ca,unsigned char**b,int rb,int cb,unsigned char*** out);
  bool GetKarnel(unsigned char** a,int r,int c,Element***out,int& rowOut);
  bool BinaryMatrixRowEchlon(unsigned char**a,int r,int col,unsigned char*** out);
  int GetCosetSize(long long int s);
  int GetHammingWt(long long int a);
  long long int GetCosetLeader(long long int a);
  public:
  int _prime;
  size_t  _len;
  unsigned char* _ir;
  Element* _one;
  Element* _zero;
  Element* _a;//primitive element
  //
  Element*** _mtrx;
  Element*** _dualBasis;
  Element*** _invMatrix;
  unsigned char** _bMatrix;
  unsigned char** _bInv;
public:
  bool Reduce(Element** e);
  bool Pad(Element** a,int l);
  bool Unpad(Element** el);
};


// c/d Tor section
Field::Field(int p,size_t l,const unsigned char* ir)
{

cout<<"\nloading...\n";
  _prime = p;
  _len = l;
  _ir = new unsigned char[_len];
  unsigned char* v = new unsigned char [_len-1];
  unsigned char* u = new unsigned char[_len-1];
  for(int i = 1; i< _len-1;i++)
    {
      v[i] = 0;
      u[i] = 0;
    }
  v[0] = 1;
  u[0] = 0;
  if(_one != 0)
    {
      delete _one;
      _one = 0;
    }
  _one = new Element(v,_len-1);
  if(_zero != 0)
    {
      delete _zero;
      _zero = 0;
    }
  _zero = new Element(u,_len-1);
  if(_a != 0)
    {
      delete _a;
      _a=0;
    }
  _a = new Element(_zero);
  _a->_val[1]= 1;
  memcpy(_ir,ir,l);
  delete []v;
  v = 0;
  delete [] u;
  u = 0;
  //
  LoadBasis();
cout<<"\n loading done..\n";

}


Field::~Field()
{
  if(_ir!= 0)
    {
      delete [] _ir;
      _ir = 0;
    }

  if(_one != 0)
    {
      delete _one;
      _one = 0;
    }

  if(_zero != 0)
    {
      delete _zero;
      _zero = 0;
    }

  if(_a != 0)
    {
      delete _a;
      _a = 0;
    }
  FreeBasis();

}

bool Field::LoadBasis()
{
int fl = _len-1;
  //Field* f = new Field(2,fl+1,(unsigned char*) "1101");// ir = 1+x+x^3
  //Element* a = new Element((unsigned char*)"010",fl);
  //
  //Forming the matrix.....
  int r = fl;
  _mtrx  =0;// new Element**[r];
  MatrixAllocate(&_mtrx,r,r);
 //
  long long int p =1;
  for(int i = 0; i< r;i++)
    {
      Element*c = 0;
      Power(_a,p,&c);
      p *= 2;
      p = (long long int) p %  (long long int)(pow(2,_len-1)-1);
      for(int j = 0 ; j < r;j++)
	{
	  Power(c,j,&_mtrx[j][i]);
	}
      delete c;
      c=0;
    }
 
  // print the matrix
  // cout<<"\n the alpha matrix\n";  
   //MatrixPrint(_mtrx,r,r);
   //cout<<"\n\n";
  //////  invert the matrix////////////////////
  _invMatrix = 0;
  MatrixInvert(_mtrx,r,&_invMatrix);
  //
  //Element*** chk = 0;
  int rr;
  int cc;
  //r = 3;
  //cout<<"i am here";
  //MatrixMultiply(_mtrx,r,r,_invMatrix,r,r,&chk,rr,cc);
  //MatrixPrint(chk,rr,cc);
  /////////////////////////////////////////////
  _dualBasis=0;
  MatrixAllocate(&_dualBasis,r,1);
  //
  for(int i = 0; i<r;i++)
    {
      _dualBasis[i][0] = new Element(_invMatrix[0][i]);
    }
  //
  cout<<"\n\ndual basis:\n\n";
  MatrixPrint(_dualBasis,r,1);
  //
  _bMatrix = 0;
  ConvertToLowerFieldMatrix(_dualBasis,r,1,&_bMatrix);
  //
  _bInv = 0;
  BinaryMatrixInvert(_bMatrix,r,&_bInv);
  //
  //cout<<"\n\nbinary matrix\n\n";
  for(int i = 0; i< r; i++)
    {
      for(int j = 0;j< r ;j++)
	{
	  int bv = _bInv[i][j];
	  //cout<<bv<<"    ";
	}
      //cout<<endl;
    }

}


bool Field::FreeBasis()
{
  //TO DO..
  int r = _len-1;
  MatrixDelete(&_mtrx,r,r);
  MatrixDelete(&_invMatrix,r,r);
  //MatrixDelete(&chk,rr,cc);
  MatrixDelete(&_dualBasis,r,1);
  //MatrixDelete(&traceCoff,1,fl);
  //
  for(int i = 0;i<r;i++)
    {
      delete [] _bMatrix[i];
      _bMatrix[i] = 0;
      //
      delete [] _bInv[i];
      _bInv[i] = 0;
      //
    }
    
  // delete []binTraceMatrix[0];
  // binTraceMatrix[0] = 0;
  // delete [] binTraceMatrix;
  // binTraceMatrix = 0;
  //
  delete [] _bMatrix;
  _bMatrix = 0;
  delete [] _bInv;
  _bInv = 0;
 
}

///!Function Unpad////////////////////////////////////////////////////////////////
// Type         : protected to Field
// Description  : Adds two elements of F[x]
// input param  : a ; b
// output param : out ...
bool Field::Unpad(Element** el)
{
  bool bRetVal = false;
  int i = 0;
  
  if((*el)->_len == 1)
    {
      return true;
    }
  for(i = (*el)->_len -1;i>=0;i--)
    {
      int vl = (*el)->_val[i];
      if(vl >= 48)
	{
	  vl = vl - 48;
	}
      if(vl == 1)
	{
	  break;
	}
    }
   
  if(i<0)
    {
      i++;
    }
  unsigned char* v = new unsigned char[i+1];
  
  memcpy(v,(*el)->_val,i+1);
  delete (*el);
  *el = 0;
  *el = new Element(v,i+1);
  delete  [] v;v=0; 
  //
  bRetVal = true;
  return bRetVal;
}

///!Function Pad////////////////////////////////////////////////////////////////
// Type         : protected to Field
// Description  : Adds two elements of F[x]
// input param  : a ; b
// output param : out ...

bool Field::Pad(Element** a,int l)
{
  bool bRetVal = false;
  if(l == 0)
    {
      return false;
    }
  unsigned char* v = new unsigned char [l];
  for(int i = 0;i<l;i++)
    {
      v[i] = (i<(*a)->_len)?(*a)->_val[i]:0;
    }

  delete (*a);
  *a = 0;
  *a = new Element(v,l);
 delete [] v;v=0; 
  bRetVal = true;
  return bRetVal;
}


///!Function PolynomialAdd////////////////////////////////////////////////////////////////
// Type         : protected to Field
// Description  : Adds two elements of F[x]
// input param  : a ; b
// output param : out ...

bool Field::PolynomialAdd(Element* a1,Element* b1, Element** out)
{
  bool bRetVal = false;
  //
  Element* a = new Element(a1);
  Element*b = new Element(b1);
  if(*out != 0)
    {
      delete (*out);
      *out = 0;
    } 

  if(a->_len < b->_len)
    {
      Pad(&a,b->_len);
    }
  else{
    if(a->_len > b->_len)
      {
	Pad(&b,a->_len);
      }
  }
  ///////////////////////////////////////////////////////
  unsigned char* temp = new unsigned char[a->_len];
  memcpy(temp,a->_val,a->_len);
  for(int i = 0;i< a->_len ; i++)
    {
      temp[i] = (((~temp[i]) & b->_val[i]) | (temp[i] & (~(b->_val[i])))) & '1';            
    }
  //
  if(*out == 0)
    {
      *out = new Element(temp,a->_len);
    }
  ////////////////////////////////////////
  if(temp != 0)
    {
      delete [] temp;
      temp = 0;
    }
  if(a!= 0)
    {
      delete a;
      a = 0;
    }
  if(b!= 0)
    {
      delete b;
      b = 0;
    }
  ////////////////////////////////////////
  bRetVal = true;  
  return bRetVal;
}

///!Function Add////////////////////////////////////////////////////////////////
// Type         : protected to Field
// Description  : Adds two elements of the finite field
// input param  : a ; b
// output param : out ...

bool Field::Add(Element* a,Element* b, Element** out)
{
  bool bRetVal = false;
  if(a->_len != b->_len)
    {
      cout<<"\n\n!Not conformable for addition\n";
      return bRetVal;
    }
  unsigned char* temp = new unsigned char[a->_len];
  memcpy(temp,a->_val,a->_len);
  for(int i = 0;i< a->_len ; i++)
    {
      temp[i] = (((~temp[i]) & b->_val[i]) | (temp[i] & (~(b->_val[i])))) & 1;            
    }
  //
  if(*out == 0)
    {
      *out = new Element(temp,a->_len);
    }
  ////////////////////////////////////////
  if(temp != 0)
    {
      delete [] temp;
      temp = 0;
    }

  ////////////////////////////////////////
  bRetVal = true;  
  return bRetVal;
}


///!Function Reduce////////////////////////////////////////////////////////////////
// Type         : private to Field
// Description  : reduces the element of input param by irreducable polynomial 
//                of the Field
// input param  : el
// output param : el ...

bool Field::Reduce(Element** el)
{
  bool bRetVal = true;
  if((*el)->_len <= _len)
    {
      //nothing to reduce
      return bRetVal;
    }
  //
  unsigned char* temp = new unsigned char[(*el)->_len - _len +1 ];
  //size_t sz = (*el)->_len;
  unsigned char* tmpEle = new unsigned char[(*el)->_len];
  memcpy(tmpEle,(*el)->_val,(*el)->_len); 
  //
  for(int k = (*el)->_len - _len ; k >= 0 ; k--)
    {
      temp[k] = tmpEle[_len-1+k];
      for(int j = _len+k-2 ; j >= k ; j--)
	{
          unsigned char foo = temp[k] & _ir[j-k];
	  tmpEle[j] =  (~tmpEle[j] & foo)| (tmpEle[j] & ~foo);
          tmpEle[j] &= 1; 
	}
    }  
  //
  delete [] (*el)->_val;
  (*el)->_val = 0;
  (*el)->_len = _len-1;
  (*el)->_val = new unsigned char[_len-1];
  memcpy((*el)->_val,tmpEle,_len-1);  
  /////////////////////////////////////////
  delete [] temp;
  temp = 0;
  delete [] tmpEle;
  tmpEle = 0;   
  return bRetVal;
}

///!Function polynomial divide////////////////////////////////////////////////////////////////
// Type         : private to Field
// Description  : reduces the element of input param by irreducable polynomial 
//                of the Field
// input param  : el
// output param : el ...

bool Field::PolynomialDivide(const Element* el ,const Element*v,Element** q,Element** r)
{
  if(v->_len == 1 )//&& v->_val[0] == 1)
    {
      *q = new Element((Element*)el);
      *r = new Element((unsigned char*)"0",1);
      return true;
    }
  bool bRetVal = true;
  if(*q != 0)
    {
      delete (*q);
      *q = 0;
    }

  if(*r != 0)
    {
      delete (*r);
      *r = 0;
    }

  //
  if(el->_len < v->_len)
    {
      //nothing to reduce
      return bRetVal;
    }
  //
  unsigned char* temp = new unsigned char[el->_len - v->_len +1 ];
  //size_t sz = (*el)->_len;
  unsigned char* tmpEle = new unsigned char[el->_len];
  memcpy(tmpEle,el->_val,el->_len); 
  //
  for(int k = el->_len - v ->_len ; k >= 0 ; k--)
    {
      temp[k] = tmpEle[v->_len-1+k];
      for(int j = v->_len+k-2 ; j >= k ; j--)
	{
          unsigned char foo = temp[k] & v->_val[j-k];
	  tmpEle[j] =  (~tmpEle[j] & foo)| (tmpEle[j] & ~foo);
          tmpEle[j] &= 1; 
	}
    }  
  
  unsigned char* rval  = new unsigned char[v->_len -1];
  memcpy(rval,tmpEle,v->_len -1);
  *r = new Element(rval,v->_len -1);
  *q = new Element(temp,el->_len - v->_len +1);
  Unpad(r);
  Unpad(q); 
 //memcpy((*el)->_val,tmpEle,_len-1);  
  /////////////////////////////////////////
  delete [] rval;
  rval = 0;
  delete [] temp;
  temp = 0;
  delete [] tmpEle;
  tmpEle = 0;   
  return bRetVal;
}



///!Function Multiply ////////////////////////////////////////////////////////////////
// Type         : private to Field
// Description  : multiply the elements of input param modulo irreducable polynomial 
//                of the Field
// input param  : a;b
// output param : out ...

bool Field::Multiply(Element* a, Element* b ,Element** out,bool domainField = true)
{
  bool bRetVal = true;
  int sz = a->_len-1 +  b->_len - 1 +1;
  if(*out != 0)
    {
      delete (*out);
      *out = 0;
    }
  *out = new Element();
  (*out)->_len = sz;
  (*out)->_val = new unsigned char[sz];
  for(int i = 0;i< sz;i++)
    {  
      (*out)->_val[i] = 0;
      for(int deg1 = 0 ; deg1 < a->_len ; deg1++)
	{
	  for(int deg2 = 0;deg2 < b->_len ; deg2++)
	    {
	      if(deg1 + deg2 != i)
		{
		  continue;
		}
              unsigned char t  =  (a->_val[deg1] & b->_val[deg2])  & 1;
	      (*out)->_val[i]  = (((~(*out)->_val[i]) & t) | ((*out)->_val[i] & ~t));
	                         
	      (*out)->_val[i] &= 1;

	    }  
	}
    }
    
  //
  if(domainField)
    {
      Reduce(out);
    }
  else{
    Unpad(out);  
  }
  return bRetVal;
}


///!Function Power ////////////////////////////////////////////////////////////////
// Type         : private to Field
// Description  : multiply the element el of input param modulo irreducable 
//                polynomial upto 2nd param p times 
// input param  : el;p
// output param : out ...
bool Field::Power(Element* el1, long long int p ,Element** out)
{
  bool bRetVal = true;
  if(p == 0)
    {
      //
      *out = new Element(_one);
      return true;
    }
  Element* el = new Element(el1);
  if(*out != 0)
    {
      delete (*out);
      *out = 0;  
    }
  *out = new Element(_one);
  Element *tmpEl = new Element(el->_val,el->_len);
  //
  while(p > 0)
    {
      if(p%2 == 1)
	{ 
          Element* tmp =0; 
	  Multiply(*out,tmpEl,&tmp);
	  delete (*out);*out = 0;
	  *out =  new Element(tmp->_val,tmp->_len);
          delete tmp;tmp = 0;
	}
      //
      Element* tmpEl2 = new Element(tmpEl->_val,tmpEl->_len);
      Element* tmp = 0;
      Multiply(tmpEl,tmpEl2,&tmp);
      delete tmpEl;tmpEl = 0;
      delete tmpEl2;tmpEl2 = 0;
      tmpEl = new Element(tmp->_val,tmp->_len);
      delete tmp;
      tmp = 0; 
      p = p/2;
    } 

  if(el!= 0)
    {
      delete el;
      el = 0;
    }
  return bRetVal;
}


///!Function Equals////////////////////////////////////////////////////////////////
// Type         : public to Field
// Description  : to check if 2 elements are equal 
// input param  : el;p
// output param : out ...

bool Field::Equals(Element*a1 , Element*b1)
{
  bool bRetVal = true;
  Element*a = new Element(a1);
  Element*b = new Element(b1);
  
  if(b->_len > a->_len)
    { 
      Pad(&a,b->_len);
    }
  if(a->_len > b->_len)
    {
       Pad(&b,a->_len);
    }
      
  for(int i = 0; i< a->_len; i++)
    {
      //
      int av = (int)a->_val[i];
      int bv = (int)b->_val[i];
      //cout <<"(\n"<<av<<" , "<< bv<< ") \n";
      if(av != bv)
	{
	  bRetVal = false;
	  break;
	}
    }

  delete a;
  a=0;
  delete b;
  b=0;
  //
  return bRetVal;
}

///!Function Invert ////////////////////////////////////////////////////////////////
// Type         : public to Field
// Description  : to find the inverse of an element 
// input param  : el;p
// output param : out ...

bool Field::Invert(Element* el,Element** out)
 {
   bool bRetVal = false;
   Element* b = new Element(el);
   Element* a = new Element(_ir,_len);
   Element* zero = new Element(_zero);
   
   if(Equals(b,zero))
     {
       *out = 0;//new Element((unsigned char*)"10000000",_len-1);
       return true;
     }
   //
 Element*one = new Element(_one);
 Element* x2 = new Element((unsigned char*)"1",1);
 Element* x1 = new Element((unsigned char*)"0",1);
 Element* y2 = new Element((unsigned char*)"0",1);
 Element* y1 = new Element((unsigned char*)"1",1);
 Element* q = 0;
 Element* r = 0;
 Element* tmp = 0;
 Element* y = 0;
 
 while(!Equals(b,zero))
   {
     if(q != 0 )
       {
	 delete q;
	 q=0;
       } 
     if(r!= 0)
       {
	 delete r;
	 r = 0;
       }
     
     Unpad(&b);
     Unpad(&a);
     PolynomialDivide((const Element*) a,(const Element*)b,&q,&r);
         
     /// x = x2 + q.x1
     Multiply(q,x1,&tmp,false);
     ///
     if(*out!= 0)
       {
	 delete (*out);
	 *out = 0;
       }
     PolynomialAdd(x2,tmp,out);
     //
     delete tmp;
     tmp = 0; 
     /// y = y2 + q.y1
     Multiply(q,y1,&tmp,false);
     PolynomialAdd(y2,tmp,&y);
     
     // 
     delete tmp;
     tmp = 0;
     //
     delete a;
     a = 0;
     a = new Element(b);
     delete b;
     b = 0;
     b = new Element(r);
     //if(Equals(b,zero))
     //{
	 //break;
     //}
     //
     delete r;
     r = 0;
     //
     delete x2;
     x2 = 0;
     x2 = new Element(x1);
     //
     delete x1;
     x1 = 0;
     x1 = new Element(*out);
     //
     delete y2;
     y2 = 0;
     y2 = new Element(y1);
     //
     delete y1;
     y1 = 0;
     y1 = new Element(y);
     //
     delete y;
     y =0;
   }

 Multiply(y2,one,out);
 ////////////////////////////////////////////////////////////////
 if(zero != 0)
   {
     delete zero;
     zero = 0;
   } 
 if(a != 0)
   {
     delete a;
     a = 0;
   }
 if(b != 0)
   {
     delete b;
     b = 0;
   }
 //
 if(x1 != 0)
   {
     delete x1;
     x1 = 0;
   }
 //
 if(x2 != 0)
   {
     delete x2;
     x2 = 0;
   }
 //
 if(y1 != 0)
   {
     delete y1;
     y1 = 0;
   }
 //
 if(y2 != 0)
   {
     delete y2;
     y2 = 0;
   }
 //
 if(q != 0)
   {
     delete q;
     q = 0;
   }
 //
 if(r != 0)
   {
     delete r;
     r = 0;
   }
 //
 if(y != 0)
   {
     delete y;
     y = 0;
   }
 //
 if(tmp != 0)
   {
     delete tmp;
     tmp = 0;
   }
 //
 if(one != 0)
   {
     delete one;
     one = 0;
   }
 //
 bRetVal = true;
 return bRetVal;
 }

///!Function MatrixInvert ////////////////////////////////////////////////////////////////
// Type         : public to Field
// Description  : to find the inverse of a matrix over field 
// input param  : a,r
// output param : out ...

bool Field::MatrixInvert(Element*** a,int r,Element**** out)
{
  bool bRetVal = false;
  Element* one = new Element(_one);
  Element* zro = new Element(_zero);
  //Element* zr = new Element((unsigned char*)"0",1);
  //Multiply(zr,_one,&zro);
  //delete zr;
  //zr = 0;
  //cout<<"\ni am zero ==>";zro->Show();
  Element*** b = 0;
  MatrixAllocate(&b,r,2*r);
  //
  for(int i = 0; i< r;i++)
    {
      for(int j = 0 ; j < r;j++)
	{
	  b[i][j] = new Element(a[i][j]);
	}
      //
      for(int k = r;k<2*r;k++)
	{                                                                       
	  if(i == k-r)
	    {
	      b[i][k] = new Element(one);
	    }
	  else
	    {
	      b[i][k] = new Element(zro);
	    }
	}
    }
  /////
  //TO DO.. convertion to rowechloni
int pivot = -1;
  for(int i = 0; i < r ; i++)
    {
      //find pivot
      pivot++;
      // get inverse
      Element * inv = 0;
///////////////////////////////////////////////////////////////
	if(Equals(b[i][pivot],_zero))
	{	
//
		int swpR = -1;
		for( int rowNo = i+1;rowNo < r;rowNo++)
		{
			if(!Equals(b[rowNo][pivot],_zero))
			{
				swpR = rowNo;
				break;
			}  
		} 
		if(swpR == -1)
		{
			if(i == r-1)
			{
				continue;
			}
			cout<<"\n singular matrix for sure \n";
			cout<<"\nrow = "<<i<<endl;
			scanf("%d",&swpR);//getch();
			return false;
		}	
		//
		for(int colNo = 0;colNo < 2*r;colNo)
		{
			Element* tmpEl =0;
			tmpEl = new Element(b[i][colNo]);
			delete b[i][colNo];
			b[i][colNo]=0;
			b[i][colNo]=new Element(b[swpR][colNo]);
			delete b[swpR][colNo];
			b[swpR][colNo] = 0;
			b[swpR][colNo] = new Element(tmpEl);
			delete tmpEl;
			tmpEl = 0;
		}

	}
//////////////////////////////////////////////////////////////
      	Invert(b[i][pivot],&inv);
      //int pivot = k;
      //multiply ith row with inverse
      for(int k = 0;k<2*r;k++)
	{
	  Element* tmp = 0;
	  Multiply(b[i][k],inv,&tmp);
	  delete b[i][k];
	  b[i][k] = 0;
	  b[i][k] = new Element(tmp);
	  delete tmp;
	  tmp = 0;
	}
      if(inv != 0)
	{
	  delete inv;
	  inv = 0;
	}
      // adjust other rows

      for(int rnew = 0; rnew < r;rnew++)
	{
	  if(rnew == i)
	    {
	      continue;
	    }
	  //
	  Element* m = new Element(b[rnew][pivot]);
	  //
          for(int l = 0;l<2*r;l++)
	    {
	      Element* tmp = 0;
	      Multiply(b[i][l],m,&tmp);
	      Element* tmp1 = 0;
	      Add(tmp,b[rnew][l],&tmp1);
	      delete b[rnew][l];
	      b[rnew][l] = 0;
	      b[rnew][l]= new Element(tmp1);
	      delete tmp;
	      tmp = 0;
	      delete tmp1;
	      tmp1 = 0; 
	    }
	  if(m != 0)
	    {
	      delete m;
	      m = 0;
	    }
	}
    }  

 //  assign the value in the output matrix
  MatrixAllocate(out,r,r);
  for(int i = 0; i< r;i++)
    {
      for(int j = 0 ; j < r;j++)
	{
	  (*out)[i][j] = new Element(b[i][j+r]);
	}
    }
  //////////////////////////////////////////////////////////////
  if(one != 0)
    {
      delete one;
      one = 0;
    }
  if(zro != 0)
    {
      delete zro;
      zro = 0;
    }
 // delete the matrix
  MatrixDelete(&b,r,2*r);
  //
  bRetVal = true;
  return bRetVal;
}

///!Function MatrixDelete ////////////////////////////////////////////////////////////////
// Type         : public to Field
// Description  : to deallocate the matrix 
// input param  : b
// output param : ...

bool Field::MatrixDelete(Element****b,int r,int c)
{
  bool bRetVal = false;
  if((*b) != 0)
    {
      for(int i = 0;i<r;i++)
	{
	  for(int j = 0; j< c;j++)
	    {
	      delete (*b)[i][j];
	      (*b)[i][j] = 0;
	    }
	  delete [] (*b)[i];
	  (*b)[i] = 0;
	}
      //
      delete [] (*b);
      *b = 0;
    }

  bRetVal = true;
  return bRetVal;
}

///!Function MatrixPrint ////////////////////////////////////////////////////////////////
// Type         : public to Field
// Description  : to print the matrix 
// input param  : mtrx
// output param : ...

bool Field::MatrixPrint(Element*** mtrx,int r,int c)
{
  bool bRetVal = false;
  cout<<endl;
  for(int i = 0; i< r;i++)
    {
      for(int j = 0; j < c ; j++)
	{
	  mtrx[i][j]->Show();cout << "       ";
	}
      cout<<endl;
    }
 cout<<"\n\n";
 bRetVal = true;
 return bRetVal;
}

///!Function MatrixAllocate ////////////////////////////////////////////////////////////////
// Type         : public to Field
// Description  : to allocate the matrix 
// input param  : out
// output param : out...

bool Field::MatrixAllocate(Element**** out,int r,int c)
{
  bool bRetVal = false;
  if(*out != 0)
    {
      delete [] (*out);
    }
*out = new Element**[r];
  for(int i = 0 ;  i < r;i++)
    {
      (*out)[i] = new Element*[c];
      for(int j = 0; j<c;j++)
	{
	  (*out)[i][j] = 0;
	}
    }
  //
  bRetVal = true;
  return bRetVal;
}


///!Function MatrixMultiply ////////////////////////////////////////////////////////////////
// Type         : public to Field
// Description  : to multiply two matrices 
// input param  : a,b
// output param : out...

bool Field::MatrixMultiply
(
 Element***a, int ar,int ac,
 Element***b,  int br,int bc,
 Element****out, int&cr,int&cc
 )
{
  bool bRetVal = false;
  if(ac != br)
    {
      cout<<"\n oops Matrices are not conformable for multiplication \n\n";
      return false;
    }
  //
  cr = ar;
  cc = bc;
    
  MatrixAllocate(out,cr,cc);
  Element* zro = new Element((unsigned char*)"0",1);
  for(int r = 0;r<cr;r++)
    {
      for(int c = 0 ; c<cc;c++)
	{
	  Multiply(zro,_one,&((*out)[r][c]));
          for(int k = 0;k<ac;k++)
	    {
	      Element* tmp = 0;
	      Multiply(a[r][k],b[k][c],&tmp);
	      //
	      Element* tmp1 = 0;
	      Add((*out)[r][c],tmp,&tmp1);
	      delete (*out)[r][c];
	      (*out)[r][c]=0;
	      (*out)[r][c] = new Element(tmp1);
	      delete tmp;tmp=0;
	      delete tmp1;tmp1=0;
	    }
	}
    }

  if(zro != 0)
    {
      delete zro;
      zro = 0;
    }
  bRetVal = true;
  return bRetVal;
}


///!Function  Matrix Power ////////////////////////////////////////////////////////////////
// Type         : public to Field
// Description  : to find a^p field element 
// input param  : a,p
// output param : out...

 bool Field::MatrixPower(Element***a,int ar,int ac,int p,Element****out, int& outr, int& outc)
  {

      bool bRetVal = false;
      ////////////////////////
      MatrixAllocate(out,ar,ac);
      outr = ar;outc=ac;
      for(int i = 0; i< ar;i++)
      {
          for(int j = 0 ; j < ac; j++)
          {
              if( i == j)
                  (*out)[i][j] = new Element(_one);
              else
                  (*out)[i][j] = new Element(_zero);
          }
      }

      if(p == 0)
      {
        //
        return true;
      }

      Element*** el = 0;
      Element***tmpEl = 0;
      MatrixAllocate(&el,ar,ac);
      MatrixAllocate(&tmpEl,ar,ac);
      for(int i = 0; i< ar;i++)
      {
            for(int j = 0 ; j < ac; j++)
            {
                el[i][j] = new Element(a[i][j]);
                tmpEl[i][j] = new Element(a[i][j]);
            }
      }
  
      int tmpr,tmpc;
      //
      while(p > 0)
      {
          /////////////////////////////////////////////////////////////////
          if(p%2 == 1)
	      {
              Element*** tmp =0; 
	          MatrixMultiply(*out,outr,outc,tmpEl,ar,ac,&tmp,tmpr,tmpc);
              for(int i = 0; i< ar;i++)
              {
                  for(int j = 0 ; j < ac; j++)
                  {
                      delete (*out)[i][j];(*out)[i][j] = 0;
                      (*out)[i][j] = new Element(tmp[i][j]);
                  }
              }
              //
              MatrixDelete(&tmp,tmpr,tmpc);
	      }
          ////////////////////////////////////////////////////////////////    
          Element*** tmpEl2=0;
          MatrixAllocate(&tmpEl2,ar,ac);
          for(int i = 0; i< ar;i++)
          {
              for(int j = 0 ; j < ac; j++)
              {
                    tmpEl2[i][j] = new Element(tmpEl[i][j]);
              }
          }
      
          Element*** tmp = 0;int tmpr,tmpc;
          MatrixMultiply(tmpEl,ar,ac,tmpEl2,ar,ac,&tmp,tmpr,tmpc);
          //MatrixDelete(&tmpEl,ar,ac);
          MatrixDelete(&tmpEl2,ar,ac);
          for(int i = 0 ; i < ar;i++)
          {
                for(int j = 0 ; j < ac;j++)
                {
                    delete tmpEl[i][j];
                    tmpEl[i][j]=0;
                    tmpEl[i][j]=new Element(tmp[i][j]);
                }
          }
          MatrixDelete(&tmp,ar,ac);
          p = p/2;
      } 

      if(el!= 0)
      {
          MatrixDelete(&el,ar,ac);
      }
      //
      if(tmpEl != 0 )
      {
          MatrixDelete(&tmpEl,ar,ac);
      }
      ////////////////////////
      bRetVal = true;
      return bRetVal;
  }



///!Function Trace ////////////////////////////////////////////////////////////////
// Type         : public to Field
// Description  : to find trace of a field element 
// input param  : a
// output param : b...

bool Field::Trace(Element*a,Element** b,int m = 1)
{
  bool bRetVal = false;
  ///////////////////////////////////////
  Element* zr = new Element((unsigned char*)"0",1);
  Multiply(zr,_one,b);
  delete zr;
  zr = 0;
  //
  int range = _len - 1-1;
  if( m > 1)
    {
      range = ((_len-1)/m) -1;
    }
  //
  for(int i = 0;i<=range;i++)
    {
      Element* tmp = 0;
      int p = pow(2,m*i);
      Power(a,p,&tmp);
      //
      Element* tmp1 = 0;
      Add((*b),tmp,&tmp1);
      delete (*b);
      *b = 0;
      *b = new Element(tmp1);
      delete tmp;tmp=0;
      delete tmp1;tmp1=0;
    }

  //////////////////////////////////////
  bRetVal = true;
  return bRetVal;
}

///!Function ConvertToLowerFieldMatrix ////////////////////////////////////////////
// Type         : public to Field
// Description  : to find trace of a field element 
// input param  : a
// output param : b...

bool Field::ConvertToLowerFieldMatrix(Element***el,int r,int c,unsigned char*** out)
{
  int row;
  int column;
  bool rowVector = true;
  if(r>1 && c >1)
    {
      *out = 0;
      return false;
    }

  if(r==1)
    {
      row = _len-1;
      column = c;
    }
  else{
    row = r;
    column = _len - 1;
    rowVector = false;
  }

  //
  *out = new unsigned char*[row];
  for(int i = 0;i<row;i++)
    {
      (*out)[i] = new unsigned char[column];
      //
      for(int j = 0; j< column;j++)
	{
	  if(rowVector)
	    {
	      (*out)[i][j] = el[0][j]->_val[i];
	    }
	  else{
	    //
	    (*out)[i][j] = el[i][0]->_val[j];
	  }
	}
    }
  //
  return true;
}


///!Function ConvertToLowerFieldMatrix ////////////////////////////////////////////
// Type         : public to Field
// Description  : to find trace of a field element 
// input param  : a
// output param : b...

bool Field::ConvertToLowerFieldMatrix(Element*el,bool rMatrix ,unsigned char*** out)
{
  int row;
  int column;
  bool rowVector = false;
  if(el == 0)
    {
      *out = 0;
      return false;
    }

  if(rMatrix == true) // its a row  matrix
    {
      row = 1;
      column = _len - 1;
      rowVector = true;
    }
  else{
    //column matrix
    row = _len - 1;
    column = 1;
      
  }
  //el->Show();
  //cout<<endl<<row<<"---------------"<<column<<endl;
  //
  *out = new unsigned char*[row];
  (*out)[0] = 0;


  for(int i = 0;i<row;i++)
    {
      (*out)[i] = new unsigned char[column];
      //
      //
      for(int j = 0; j< column;j++)
	{
	  if(rowVector)
	    {
	      //
	      //cout<<"\nasssssssssssssssssssssinnn\n";el->Show();
	      (*out)[i][j] = el->_val[j];
	      //cout<<"\n doneeeeeeeeeeeeeee\n\n";
	    }
	  else{
	    //
	    (*out)[i][j] = el->_val[i];
	  }
	}
    }
  //

  return true;
}




bool Field::BinaryMatrixInvert(unsigned char**a,int r,unsigned char*** out)
{
  unsigned char** temp = new unsigned char*[r];
  for(int i = 0; i< r;i++)
    {
      temp[i] = new unsigned char[2*r];
    }
  //
  //
  for(int i = 0; i<r;i++)
    {
      for(int j = 0; j< r;j++)
	{
	  temp[i][j] = a[i][j];
	}
      //
      for(int k = r;k<2*r;k++)
	{
	  if(k-r == i)
	    {
	      temp[i][k] = 1;
	    }
	  else{
	    temp[i][k] = 0;
	  }
	}
    }
  ////////////////////////////////////////////////
  // convert into row echlon
  int pivot = -1;
  for(int i = 0; i< r;i++)
    {
      pivot++;
      
      // swap rows if needed//////
      if(temp[i][pivot] == 0)
	{ int swpR = -1;
	  for(int s = i;s<r;s++)
	    {
	      if(temp[s][pivot] !=0)
		{
		 swpR = s;
		 break;
		}
	    }
	    if(swpR == -1)
		{
                  //continue;
		  int swpCol = -1;
                  for(int colWise = pivot;colWise < r;colWise++)
		    {
		      if(temp[i][colWise] != 0)
			{
			  //swpR = colWise;
			  pivot = colWise;
			  break;
			}
		    }
		  if(swpCol == -1)
		    {
		      cout<< "\nIts a singular matrix\n\n";
		      // put the zero row at  the bottom
		      return false;
		    }
		}
	      //
	    if(swpR != -1)
	      {
	      for(int k = 0;k< 2*r;k++)
		{
		  unsigned char t =temp[i][k];
		  temp[i][k] = temp[swpR][k];
		  temp[swpR][k]=t;
		}
	      }
	 }/////////////////////////////////////////////
	  //
	  for(int k = 0;k<r;k++)
	    {
	      if(k == i)
		{
		  continue;
		}
	      //
	      if(temp[k][pivot] == 1)
		{
		  for(int l = 0; l< 2*r;l++)
		    {
		      temp[k][l] = (~temp[k][l] & temp[i][l]) |  (temp[k][l] & ~temp[i][l]);
		      temp[k][l] &= 1;  
		    }
		}
	    }
	}

      /////////////////////////////
  *out = new unsigned char*[r];
  for(int i = 0; i < r;i++)
    {
      (*out)[i] = new unsigned char[r];
      //
      for(int k = 0;k<r;k++)
	{
	  (*out)[i][k] = temp[i][k+r];
	}
    }
  

// To delete the temp array.......
  for(int i = 0 ; i < r ; i++)
    {
      delete [] temp[i];
      temp[i] = 0;
    }
  delete [] temp;
}


bool Field::BinaryMatrixRowEchlon(unsigned char**a,int r,int col,unsigned char*** out)
{
  unsigned char** temp = new unsigned char*[r];
  for(int i = 0; i< r;i++)
    {
      //
      temp[i] = 0;
      temp[i] = new unsigned char[col];
    }
  //
  //
  //cout<<"\ntemp:\n";
  for(int i = 0; i<r;i++)
    {
      for(int j = 0; j< col;j++)
	{
	  temp[i][j] = a[i][j];
	  //int vv = (int)temp[i][j];
	  //cout<<vv<<"    ";
	}
      //
      //cout<<endl;
    }
  ////////////////////////////////////////////////
  // convert into row echlon
  int pivot = -1;
  for(int i = 0; i< r;i++)
    {
      pivot++;
      
      // swap rows if needed//////
      if(temp[i][pivot] == 0)
	{ int swpR = -1;
	  for(int s = i+1;s<r;s++)
	    {
	      if(temp[s][pivot] !=0)
		{
		 swpR = s;
		 break;
		}
	    }
	    if(swpR == -1)
		{
                   continue;
		  // int swpCol = -1;
                  // for(int colWise = pivot;colWise < col;colWise++)
		  //   {
		  //     if(temp[i][colWise] != 0)
		  // 	{
		  // 	  //swpR = i;
		  // 	  swpCol = colWise;
		  // 	  pivot = colWise;
		  // 	  break;
		  // 	}
		  //   }
		  // if(swpCol == -1)
		  //   {
		  //     cout<< "\nIts a singular matrix\n\n";
		  //     continue;
		      // //
		      // for(int rr = i+1;rr < r;rr++)
		      // 	{
		      // 	  for(int cc = 0; cc < col;cc++)
		      // 	    {
		      // 	      if(temp[rr][cc] !=0)
		      // 		{
		      // 		  swpR = rr;
		      // 		  break;
		      // 		}
		      // 	    }
		      // 	  if(swpR != -1)
		      // 	    {
		      // 	      break;
		      // 	    }
		      // 	}
		      
		   // }////////////////////////////////////////////////////
		}
	    //
	    // if(swpR == -1)
	    //   {
	    // 	continue;
            //     cout<<"\nnullllllllllllllllllll\n";
	    // 	//return false;
	    //   }
	    if(swpR != i)
	      {
		for(int k = 0;k< col;k++)
		  {
		    unsigned char t =temp[i][k];
		    temp[i][k] = temp[swpR][k];
		    temp[swpR][k]=t;
		  }
	      }
	 }/////////////////////////////////////////////
	  //
	  for(int k = 0;k<r;k++)
	    {
	      if(k == i)
		{
		  continue;
		}
	      //
	      if(temp[k][pivot] == 1)
		{
		  for(int l = 0; l< col;l++)
		    {
		      temp[k][l] = (~temp[k][l] & temp[i][l]) |  (temp[k][l] & ~temp[i][l]);
		      temp[k][l] &= 1;  
		    }
		}
	    }
    }

  /////////////////////////////
  *out = new unsigned char*[r];
  for(int i = 0; i < r;i++)
    {
      (*out)[i] = new unsigned char[col];
      //
      for(int k = 0;k<col;k++)
	{
	  (*out)[i][k] = temp[i][k];
	}
    }
  

// To delete the temp array.......
  for(int i = 0 ; i < r ; i++)
    {
      delete [] temp[i];
      temp[i] = 0;
    }
  delete [] temp;
}


bool Field::BinaryMatrixMultiply(unsigned char**a,int ra,int ca,unsigned char**b,int rb,int cb,unsigned char*** out)
{

  if(ca != rb)
    {
      cout<<"\n\nMatrices are not conformable for multiplication\n\n";
      return false;
    }
  //
  *out = new unsigned char*[ra];
  for(int i = 0; i < ra;i++)
    {
      (*out)[i]= new unsigned char[cb];
    }
  //
  for(int i = 0;i<ra;i++)
    {
      for(int j = 0 ; j < cb;j++)
	{
	  //
	  (*out)[i][j] = 0;
	  for(int k = 0; k< ca;k++)
	    {
	      unsigned char t = (a[i][k]&b[k][j])&1;
	      (*out)[i][j] = (*out)[i][j] & ~t | ~(*out)[i][j] & t;
	      //(*out)[i][j] &= 1;
	    }
	  //
	}
    }

  return true;
}



bool Field::GetKarnel(unsigned char** a,int row,int col,Element***out,int& rowOut)
{
  bool bRetVal = false;
  unsigned char** temp = 0;
  BinaryMatrixRowEchlon(a,row,col,&temp);
  //
  cout<<"\nrow echlon matrix\n\n";
  for(int i = 0; i<row; i++)
    {
      for(int j = 0 ; j<col;j++)
  	{
	  temp[i][j] &= 1;
  	  int vv = (int) temp[i][j];
  	  cout<<vv<<"    ";
  	}
      cout<<endl;
    }
  
  int** freeVariables = new int*[col];
  for(int i = 0; i< col;i++)
    {
      freeVariables[i] =new int[col];
      //
      for(int j = 0; j< col;j++)
	{
	  freeVariables[i][j] = -1; // free variables getting initialized by -1
	}
    }
  //
  //
  for(int r = row -1; r >= 0 ;r--)
    {
      bool alZro = true;
      int index = -1;
      //
      for(int c = col -1;c >= 0;c--)
	{
	  if(temp[r][c] != 0)
	    {
	      alZro = false;
	      index = c;
	      break;
	    }
	}
      //
      if(alZro)
	{
	  continue;
	}
      /////////////////////////////////////////
      // 
      bool oneNonZero = true;
      int variableIndex = index+1;
      freeVariables[variableIndex-1][variableIndex-1] = variableIndex; 
      for(int c = index-1;c >= 0;c--)
	{
	  if(temp[r][c] != 0)
	    {
	      oneNonZero = false;
	      freeVariables[c][c] = c+1;
	      //variableIndex =  c+1;
	    }
	}
      //
      if(oneNonZero == true)
	{
	  if(index == 6)
	    {
	      cout<<"\ndoooooooneeeee\n";
	    }
          //cout<<"\nthis is wrong\n";
	  for(int c = 0 ; c < col ; c++)
	    {
	      freeVariables[c][index] = 0;
	      freeVariables[index][c] = -1;
	    }
	  //
	  for(int c = 0; c < r;c++)
	    {
	      temp[c][index] = 0;
	    }
	  continue;
	}
      
      //////////////////////////
      /// fix the last variable
      int firstIndex = -1;
      for(int c = 0;c<col;c++)
	{
	  if(temp[r][c] == 1)
	    {
	      firstIndex = c;
	      break;
	    }
	}
      
      //
      for(int c = firstIndex+1;c < col;c++)
	{
	  if(temp[r][c] == 1)
	    {
	      freeVariables[c][firstIndex] = c+1;
	    }
	}
      //
      for(int fi = 0; fi < col;fi++)
	{
	  freeVariables[firstIndex][fi] = -1;
	}
    }// end of step1 of freeVariableMartix formation
  
  //
  for(int c =0;c<col;c++)
    {
      bool empty = true;
     
      for(int r = 0; r< col;r++)
	{
	  if(freeVariables[r][c] != -1)
	    {
	      empty = false;
	      break;
	    }
	}
      //
      if(empty)
	{
	  freeVariables[c][c] = c+1;
	}
    }// finally freevariable matrix formation ends
  //

  ////print freeVariables
  cout<<"\nfree variables \n\n";
  for(int i = 0 ;i<col;i++)
    {
      for(int j = 0 ; j < col;j++)
  	{
  	  cout<<freeVariables[i][j]<<"  ,  ";
  	}
      cout<<endl;
    }
  cout<<endl;
  // to form the basis
  // size of basis will be the number of nonempty rows in freeVariables array
  // let's count that
  int rowSize = 0;
  //
  int* nonEmptyRowIndex = new int[col];
  //
  for(int i = 0;i<col;i++)
    {
      bool empty = true;
      for(int j = 0;j<col;j++)
	{
	  if(freeVariables[i][j] > 0)
	    {
	      empty = false;
              
	      break;
	    }
	}
      if(!empty)
	{
          nonEmptyRowIndex[rowSize]= i;
	  rowSize++;
	}
    } 
  rowOut = rowSize;
  //
  //cout<<"\n---------------->"<<rowSize<<endl;
  /////////////////////////////////
  *out = new Element*[rowSize];
  for(int i = 0; i< rowSize;i++)
    {
      (*out)[i] = new Element(_zero);
    }
  //
  for(int i = 0; i < rowSize;i++)
    {
      int rowValue = nonEmptyRowIndex[i];
      for(int j = 0; j<col;j++)
	{
	  if(freeVariables[rowValue][j] > 0)
	    {
	      Element* ap =0;
	      Power(_a,j,&ap);
	      Element* bp = 0;
	      Add((*out)[i],ap,&bp);
	      delete (*out)[i];
	      (*out)[i] = 0;
	      (*out)[i]= new Element(bp);
	      delete ap;
	      ap = 0;
	      delete bp;
	      bp = 0;
	    }
	}
    }

  ///////////////////////////////////////
  for(int i = 0; i < col;i++)
    {
      delete [] freeVariables[i];
      freeVariables[i] = 0;
    }
  delete [] freeVariables;
  freeVariables = 0;
  //
  delete [] nonEmptyRowIndex;
  nonEmptyRowIndex = 0;
  //
  for(int i = 0; i < row;i++)
    {
      delete [] temp[i];
      temp[i] = 0;
    }
}

int Field::GetCosetSize(long long int s)
{
  s  = s % ((int) pow(2,_len-1)-1);
  if(s == 0)// || (s % pow(2,_len-1) - 1) == 0 )
    {
  return 1;
    }
  int n = pow(2,_len-1)-1;
  int sz = 0;
  int initialS = s;
  //int newS = 0;//s * 2 % n;
do
  {
    sz++;
    s = (s * 2) % n;
    //s = newS;
  }
 while(initialS != s);
 return sz;
}

long long int Field::GetCosetLeader(long long int a)
{
  int n = pow(2,_len-1)-1;
  int sz = 0;
  a = a % (int) (pow(2,_len-1)-1);
  int initialA = a;
  if(a == 0)
    {
      return 0;
    }
  int* ary = new int[_len-1];
  //int newS = 0;//s * 2 % n;
do
  {
    sz++;
    a = (a * 2) % n;
    ary[sz-1]= a;
  }
 while(initialA != a);
 //
  for(int i = 0; i <  sz;i++)
    {
      for(int j = sz-1;j>i;j--)
	{
	  if((ary[j]) < (ary[j-1]))
	    {
	      int tmp = ary[j-1];
	      ary[j-1]=ary[j];
	      ary[j]=tmp;
	    }
	}
    }
  //
  int v = ary[0];
  delete ary;
  ary = 0;
  return v;
}

int Field::GetHammingWt(long long int a)
{
  int wt = 0;
  while(a>0)
    {
      if(a%2 == 1)
	{
	  wt++;
	}
      a = a/2;
    }
  return wt;
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//!                         TEST BED                                      !//
/////////////////////////////////////////////////////////////////////////////
