class Matrix:
    "Generate and represent matrix objects"


    """Initialize a matrix object with at least one row and one column, and the same number of columns in each row."""
    def __init__(self, matrix):
        self.rows = len(matrix)
        if not self.rows:
            raise ValueError("Matrix must have at least one row")
        self.cols = len(matrix[0])
        if not self.cols:
            raise ValueError("Matrix must have at least one column")
        self.matrix = tuple(tuple(map(float, row)) for row in matrix if len(row) == self.cols)
        if not len(self.matrix) == self.rows:
            raise ValueError("All rows must have the same number of columns")
        self.square = self.rows == self.cols


    @staticmethod
    def sign(row, col):
        return -1 if (row + col) % 2 else 1

    def submatrix(self, row, col):
        return Matrix(tuple(tuple(self.matrix[i][j]
            for j in range(self.cols) if j != col)
            for i in range(self.rows) if i != row))

    def determinant(self):
        if not self.square:
            raise ValueError("Matrix must be square")
        if self.rows == 1:
            return self.matrix[0][0]
        return sum(Matrix.sign(i,0)*
            self.matrix[i][0]*self.submatrix(i,0).determinant()
            for i in range(self.rows))
    
    def __add__(self, other):
        if not isinstance(other, Matrix):
            raise TypeError("Can only add matrices")
        if not self.rows == other.rows or not self.cols == other.cols:
            raise ValueError("Matrices must be of the same size")
        return Matrix(tuple(tuple(self.matrix[i][j] + other.matrix[i][j]
            for j in range(self.cols))
            for i in range(self.rows)))

    __radd__ = __add__