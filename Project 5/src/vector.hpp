

//! This class holds an array of type T
/*! 
 * This class supports addition, subtraction and 
 * tridiagonal matrix operations on simple matrices
 * 
 */
template <typename T>
class Vector
{
private:
  //! The size of the vector
  int N;
  
  //! The array elements
  T* vec;
  
public:
  //! Creates the vector
  Vector(const int size) : N(size)
  {
    vec = new T[N];
  }
  
  //! Destructs the vector
  ~Vector()
  {
    delete[] vec;
  }
  
  //! Accesing of the elements
  T& operator[](const int N)
  {
    return vec[N];
  }
  
  //! Constant version of accessor
  T operator[](const int N) const
  {
    return vec[N];
  }
  
  //! Solves the equation d = Au
  /*!
   * d is this vector (known), A is on the form 
   * \f$
   * \begin{pmatrix}
   * a & b & 0 & 0\\
   * b & a & b & 0\\
   * 0 & b & a & b
   * \end{pmatrix}\f$
   * and u is the returned vector
   */
  Vector<T> solve(const Vector &d, const T a, const T b) const
  {
    Vector<T> u(N);
    
    Vector<T> beta(N);
    Vector<T> gamma(N);
    
    beta[0] = vec[0]/a;
    gamma[0] = -b/a;
    
    for (int ii=1; ii<N; ii++){
      beta[ii] = (d[ii] - b*beta[ii-1])/(b*gamma[ii-1] + a);
      gamma[ii] = -b/(b*gamma[ii-1]+a);
    }
    u[N-1] = beta[N-1];
    for (int ii=N-1; ii>1; ii--){
      u[ii-1] = beta[ii-1] + gamma[ii-1]*u[ii];
    }
    
    return u;
  }
  
  //! Solves d = Au
  /*!
   * A is known, and u is this vector. A is on the form of solve()
   * this is done destructively? See if destructive is necessary/wanted
   * 
   */
  Vector<T> multiply(const T a, const T b) const
  {
    Vector<T> d(N);
    d[0] = a*vec[0] + b*vec[1];
    d[N-1] = b*vec[N-2] + a*vec[N-1];
    
    for (int ii=1; ii<N-1; ii++){
      d[ii] = b*vec[ii-1] + a*vec[ii] + b*vec[ii+1];
    }
   return d;
  }
  
};