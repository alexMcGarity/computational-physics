class Matrix :
  "Generate and manipulate matrix objects"

  def __init__ (self,matrix) : # Object constructor
    "initialize a matrix object"
    self.rows = len( matrix )
    if not ( self.rows ) :
      raise ValueError("Invalid Matrix")
    self.cols = len( matrix[0] )
    if not ( self.cols ) :
      raise ValueError("Invalid Matrix")
    self.matrix = tuple( tuple( map( float, row ))
      for row in matrix if ( len(row) == self.cols ))
    if not ( len(self.matrix) == self.rows ) :
      raise ValueError("Invalid Matrix")
    self.square = ( self.cols == self.rows )

  @staticmethod
  def unit (size,val=1.0) :
    "returns a unit matrix of a given size"
    return Matrix( tuple( tuple(
      ( float(val) if ( i == j ) else 0.0 )  
      for j in range(size))
      for i in range(size)))

  @staticmethod
  def sign (row,col) :
    "returns the sign in the signed cofactor of a matrix"
    return ( -1.0 if ( int( row + col ) % 2 ) else +1.0 )

  def submatrix (self,row,col) :
    "returns a submatrix of a matrix"
    return Matrix( tuple( tuple( self.matrix[i][j]
      for j in range(self.cols) if ( j != col ))
      for i in range(self.rows) if ( i != row )))

  def determinant (self) :
    "returns the determinant of a matrix"
    if not ( self.square ) :
      raise ValueError("Non-Square Matrix")
    if ( self.rows == 1 ) :
      return self.matrix[0][0]
    return sum(( Matrix.sign(i,0) *
      self.matrix[i][0] *
      self.submatrix(i,0).determinant())
      for i in range(self.rows))

  def __add__ (self,other) :
    "returns the sum of two matrices or a matrix and a scalar"
    if not isinstance( other, Matrix ) :
      return ( self + Matrix.unit( self.cols, other ))
    if (( other.rows != self.rows ) or
      ( other.cols != self.cols )) :
      raise ValueError("Invalid Matrix Sum")
    return Matrix( tuple( tuple(
      ( self.matrix[i][j] + other.matrix[i][j] )
      for j in range(other.cols))
      for i in range(self.rows)))

  __radd__ = __add__

  def __getitem__ (self,key) :
    "returns a row of a matrix"
    return self.matrix[key]

  def __str__ (self) :
    "returns a string representation of a matrix"
    return "\n".join( " ".join( format( val, "+9.2E" )
      for val in row )
      for row in self.matrix )
  
  def trace (self) :
    "returns the trace of a matrix"
    if not ( self.square ) :
      raise ValueError("Non-Square Matrix")
    return sum( self.matrix[i][i]
      for i in range(self.rows))
  
  def transpose (self) :
    "returns the transpose of a matrix"
    return Matrix( tuple( tuple( self.matrix[i][j]
      for i in range(self.rows))
      for j in range(self.cols)))
  
  def __sub__ (self,other) :
    "returns the difference of two matrices or a matrix and a scalar"
    if not isinstance( other, Matrix ) :
      return ( self - Matrix.unit( self.cols, other ))
    if (( other.rows != self.rows ) or
      ( other.cols != self.cols )) :
      raise ValueError("Invalid Matrix Sum")
    return Matrix( tuple( tuple(
      ( self.matrix[i][j] - other.matrix[i][j] )
      for j in range(other.cols))
      for i in range(self.rows)))
  
  def __mul__ (self,other) :
    "returns the product of two matrices or a matrix and a scalar"
    if not isinstance( other, Matrix ) :
      return ( self * Matrix.unit( self.cols, other ))
    if ( self.cols != other.rows ) :
      raise ValueError("Invalid Matrix Product")
    return Matrix( tuple( tuple(
      sum(( self.matrix[i][k] *
        other.matrix[k][j] )
        for k in range(self.cols))
      for j in range(other.cols))
      for i in range(self.rows)))
  

  
  def inverse (self) :
    "returns the inverse of a matrix using the determinant method"
    if not ( self.square ) :
      raise ValueError("Non-Square Matrix")
    det = self.determinant()
    if ( det == 0.0 ) :
      raise ValueError("Singular Matrix")
    return ( ( 1.0 / det ) *
      Matrix( tuple( tuple(
        ( Matrix.sign(i,j) *
          self.submatrix(i,j).determinant())
          for j in range(self.cols))
          for i in range(self.rows))).transpose())