import random


class Matrix :
  "Generate and manipulate matrix objects"

  def __init__ (self,matrix) : # Object constructor
    "initialize a matrix object"
    self.rows = len( matrix )
    # Check the rows and columns of the matrix for validity
    if not ( self.rows ) :
      raise ValueError("Invalid Matrix")
    self.cols = len( matrix[0] )
    if not ( self.cols ) :
      raise ValueError("Invalid Matrix")
    self.matrix = tuple( tuple( map( float, row ))
      for row in matrix if ( len(row) == self.cols ))
    # throw an exception if the matrix is not rectangular
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
  
  def reducedRowEchelonForm (self) :
    "New Function: returns the reduced row echelon form of a matrix using Gaussian Elimination with partial pivoting, returns the reduced matrix and the number of row swaps"
    # Initialize the reduced matrix and the number of row swaps
    reducedMatrix = list( list( self.matrix[i][j]
        for j in range(self.cols))
        for i in range(self.rows))
    rowSwaps = 0
    # Iterate through the rows of the matrix
    for i in range( self.rows ) :
      # Find the maximum value in the current column
      maxValue = abs( reducedMatrix[i][i] )
      maxIndex = i
      for j in range( i + 1, self.rows ) :
        if ( abs( reducedMatrix[j][i] ) > maxValue ) :
          maxValue = abs( reducedMatrix[j][i] )
          maxIndex = j
      # Swap the current row with the row containing the maximum value
      if ( maxIndex != i ) :
        reducedMatrix[i], reducedMatrix[maxIndex] = reducedMatrix[maxIndex], reducedMatrix[i]
        rowSwaps += 1
      # Iterate through the remaining rows
      for j in range( i + 1, self.rows ) :
        # Check for a zero pivot
        if ( reducedMatrix[i][i] == 0.0 ) :
          raise ValueError("Zero Pivot")
        # Calculate the multiplier
        multiplier = reducedMatrix[j][i] / reducedMatrix[i][i]
        # Iterate through the columns
        for k in range( self.cols ) :
          reducedMatrix[j][k] -= multiplier * reducedMatrix[i][k]
    # Return the reduced matrix and the number of row swaps
    return reducedMatrix, rowSwaps

  def determinant (self) :
    "returns the determinant of a matrix using the reduced row echelon form to accelerate the calculation"
    if not ( self.square ) :
      raise ValueError("Non-Square Matrix")
    # Calculate the reduced row echelon form of the matrix
    reducedMatrix, rowSwaps = self.reducedRowEchelonForm()
    # Calculate the determinant
    determinant = 1.0
    for i in range( self.rows ) :
        determinant *= reducedMatrix[i][i]
    # Adjust the sign of the determinant based on the number of row swaps
    if ( rowSwaps % 2 ) :
        determinant *= -1.0
    # Return the determinant
    return determinant

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
    # throw an exception if the matrices are not the same size
    if (( other.rows != self.rows ) or
      ( other.cols != self.cols )) :
      raise ValueError("Invalid Matrix Sum")
    
    # The following code is equivalent to the following
    # for i in range(self.rows) :
    #   for j in range(other.cols) :
    #     result[i][j] = self.matrix[i][j] - other.matrix[i][j]
    
    return Matrix( tuple( tuple(
      ( self.matrix[i][j] - other.matrix[i][j] )
      for j in range(other.cols))
      for i in range(self.rows)))
  
  __rsub__ = __sub__
  
  def __mul__ (self,other) :
    "returns the product of two matrices or a matrix and a scalar"
    if not isinstance( other, Matrix ) :
      return ( self * Matrix.unit( self.cols, other ))
    if ( self.cols != other.rows ) :
      raise ValueError("Invalid Matrix Product")
    
    # The following code is equivalent to the following
    # for i in range(self.rows) :
    #   for j in range(other.cols) :
    #     for k in range(self.cols) :
    #       result[i][j] += self.matrix[i][k] * other.matrix[k][j]

    return Matrix( tuple( tuple(
      sum(( self.matrix[i][k] *
        other.matrix[k][j] )
        for k in range(self.cols))
      for j in range(other.cols))
      for i in range(self.rows)))
  
  __rmul__ = __mul__
  
  def inverse (self) :
    "returns the inverse of a matrix using Gaussian Elimination with partial pivoting"
    if not ( self.square ) :
      raise ValueError("Non-Square Matrix")
    if ( self.determinant() == 0.0 ) :
        raise ValueError("Singular Matrix")
    # suppose the matrix is A, then we want to find the inverse of A
    # we will use the augmented matrix [A|I] and perform row operations
    # to transform [A|I] into [I|A^-1]
    # first, create the augmented matrix
    augmentedMatrix = list( list( self.matrix[i][j]
        for j in range(self.cols))
        for i in range(self.rows))
    for i in range(self.rows) :
        augmentedMatrix[i] += list( 0.0
            for j in range(self.cols))
        augmentedMatrix[i][i+self.cols] = 1.0
    # perform row operations to transform the augmented matrix into the identity matrix
    for i in range(self.rows) :
        # find the row with the largest pivot
        maxIndex = i
        for j in range( i + 1, self.rows ) :
            if ( abs( augmentedMatrix[j][i] ) > abs( augmentedMatrix[maxIndex][i] ) ) :
                maxIndex = j
        # swap the current row with the row containing the largest pivot
        augmentedMatrix[i], augmentedMatrix[maxIndex] = augmentedMatrix[maxIndex], augmentedMatrix[i]
        # divide the current row by the pivot
        pivot = augmentedMatrix[i][i]
        for j in range( 2 * self.cols ) :
            augmentedMatrix[i][j] /= pivot
        # eliminate the pivot from the other rows
        for j in range( self.rows ) :
            if ( j != i ) :
                multiplier = augmentedMatrix[j][i]
                for k in range( 2 * self.cols ) :
                    augmentedMatrix[j][k] -= multiplier * augmentedMatrix[i][k]
    # return the inverse of the matrix
    return Matrix( tuple( tuple( augmentedMatrix[i][j]
        for j in range(self.cols, 2 * self.cols))
        for i in range(self.rows)))
  
  def solve (self, Y):
    "Solves the system of equations for the matrix and the vector Y"
    if not isinstance( Y, Matrix ) :
      raise ValueError("Invalid Matrix")
    if ( Y.rows != self.rows ) :
        raise ValueError("Invalid Matrix")
    if ( Y.cols != 1 ) :
        raise ValueError("Invalid Matrix")
    return ( self.inverse() * Y )
  
    
