#ifndef MATRIX_SRC_S21_MATRIX_OOP_H
#define MATRIX_SRC_S21_MATRIX_OOP_H

#include <iostream>

class S21Matrix {
 private:
  // Attributes
  int rows_, cols_;  // Rows and columns
  double **matrix_;  // Pointer to the memory where the matrix is allocated

  const float minimum_diff_ = 1e-7;

  [[nodiscard]] S21Matrix minor(int m, int n) const;
  [[nodiscard]] bool row_column_equal(const S21Matrix &A) const noexcept;
  [[nodiscard]] double mnoj_matrix(const S21Matrix &A, int i,
                                   int j) const noexcept;
  [[nodiscard]] double determinant_out() const;
  void CreateMatrix();
  void CopyMatrix(const S21Matrix &A);
  void DeleteMatrix(S21Matrix &A) noexcept;

 public:
  S21Matrix() noexcept;
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other) noexcept;
  ~S21Matrix();

  [[nodiscard]] inline int getRows() const noexcept;
  void setRows(int input);
  [[nodiscard]] int getCols() const noexcept;
  void setCols(int input);

  [[nodiscard]] bool EqMatrix(const S21Matrix &other) const noexcept;
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num) noexcept;
  void MulMatrix(const S21Matrix &other);
  [[nodiscard]] S21Matrix Transpose() const;
  [[nodiscard]] S21Matrix CalcComplements() const;
  [[nodiscard]] double Determinant() const;
  [[nodiscard]] S21Matrix InverseMatrix() const;

  S21Matrix operator+(const S21Matrix &A) const;
  S21Matrix operator-(const S21Matrix &A) const;
  S21Matrix operator*(const double number) const noexcept;
  S21Matrix operator*(const S21Matrix &A) const;
  bool operator==(const S21Matrix &A) const noexcept;
  bool operator!=(const S21Matrix &A) const noexcept;
  S21Matrix &operator=(S21Matrix &A);
  S21Matrix &operator=(S21Matrix &&A) noexcept;
  S21Matrix operator+=(const S21Matrix &A);
  S21Matrix operator-=(const S21Matrix &A);
  S21Matrix operator*=(const S21Matrix &A);
  S21Matrix operator*=(const double number) noexcept;
  double &operator()(int i, int j);
  const double &operator()(int i, int j) const;
  friend S21Matrix operator*(const double num, const S21Matrix &A) noexcept;
};

#endif  // MATRIX_SRC_S21_MATRIX_OOP_H