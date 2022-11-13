#include "s21_matrix_oop.h"

// Constructors

S21Matrix::S21Matrix() noexcept : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ <= 0 || cols_ <= 0)
    throw std::out_of_range(
        "Incorrect input, rows and cols size should be positive");
  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
  CopyMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix &&other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() {
  if (matrix_) {
    DeleteMatrix(*this);
  }
}

// accessors and mutators

int S21Matrix::getRows() const noexcept { return rows_; }

void S21Matrix::setRows(int input) {
  if (input <= 0)
    throw std::out_of_range("Incorrect input, size should be positive");
  if (input != rows_) {
    double **sol = new double *[input];
    for (int i = 0; i < input; i++) {
      sol[i] = new double[cols_];
      for (int j = 0; j < cols_; j++) {
        if (i < rows_) {
          sol[i][j] = matrix_[i][j];
        } else {
          sol[i][j] = 0;
        }
      }
    }
    DeleteMatrix(*this);
    rows_ = input;
    matrix_ = sol;
  }
}

int S21Matrix::getCols() const noexcept { return cols_; }

void S21Matrix::setCols(int input) {
  if (input <= 0)
    throw std::out_of_range("Incorrect input, size should be positive");
  if (input != cols_) {
    double **sol = new double *[rows_];
    for (int i = 0; i < rows_; i++) {
      sol[i] = new double[input];
      for (int j = 0; j < input; j++) {
        if (j < cols_) {
          sol[i][j] = matrix_[i][j];
        } else {
          sol[i][j] = 0;
        }
      }
    }
    DeleteMatrix(*this);
    cols_ = input;
    matrix_ = sol;
  }
}

// private functions

void S21Matrix::CreateMatrix() {
  matrix_ = new double *[rows_]();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

void S21Matrix::CopyMatrix(const S21Matrix &A) {
  rows_ = A.rows_;
  cols_ = A.cols_;
  CreateMatrix();
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = A.matrix_[i][j];
    }
  }
}

void S21Matrix::DeleteMatrix(S21Matrix &A) noexcept {
  for (int i = 0; i < A.rows_; i++) {
    delete[] A.matrix_[i];
  }
  delete[] A.matrix_;
}

S21Matrix S21Matrix::minor(int m, int n) const {
  int flagi = 0, flagj = 0;
  S21Matrix result(rows_ - 1, cols_ - 1);
  for (int i = 0; i < result.rows_; i++) {
    if (i == m) flagi = 1;
    for (int j = 0; j < result.cols_; j++) {
      if (j == n) flagj = 1;
      result(i, j) = matrix_[i + flagi][j + flagj];
    }
    flagj = 0;
  }
  return result;
}

bool S21Matrix::row_column_equal(const S21Matrix &A) const noexcept {
  return A.cols_ == cols_ && A.rows_ == rows_;
}

double S21Matrix::mnoj_matrix(const S21Matrix &A, int i, int j) const noexcept {
  double sol = 0;
  for (int k = 0; k < cols_; k++) {
    sol += matrix_[i][k] * A.matrix_[k][j];
  }
  return sol;
}

double S21Matrix::determinant_out() const {
  double result = 0;
  if (rows_ == 1) {
    result = matrix_[0][0];
  } else {
    for (int i = 0; i < rows_; i++) {
      S21Matrix buff;
      result +=
          ((i % 2) ? -1 : 1) * matrix_[0][i] * minor(0, i).determinant_out();
    }
  }
  return result;
}

// public functions

bool S21Matrix::EqMatrix(const S21Matrix &other) const noexcept {
  if (!row_column_equal(other)) return false;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (std::abs(matrix_[i][j] - other.matrix_[i][j]) >= minimum_diff_) {
        return false;
      }
    }
  }
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (!row_column_equal(other))
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (!row_column_equal(other))
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) noexcept {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] * num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) noexcept {
  if (other.rows_ != cols_)
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  S21Matrix buff(rows_, other.cols_);
  for (int i = 0; i < buff.rows_; i++) {
    for (int j = 0; j < buff.cols_; j++) {
      buff.matrix_[i][j] = mnoj_matrix(other, i, j);
    }
  }
  *this = std::move(buff);
}

S21Matrix S21Matrix::Transpose() const {
  S21Matrix sol(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      sol.matrix_[j][i] = matrix_[i][j];
    }
  }
  return sol;
}

double S21Matrix::Determinant() const {
  if (rows_ != cols_)
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  return determinant_out();
}

S21Matrix S21Matrix::CalcComplements() const {
  S21Matrix buff(rows_, cols_);
  if (rows_ != cols_)
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      buff(i, j) = (((i + j) % 2) ? -1 : 1) * minor(i, j).Determinant();
    }
  }
  return buff;
}

S21Matrix S21Matrix::InverseMatrix() const {
  double det = Determinant();
  if (std::abs(det) < minimum_diff_) throw std::out_of_range("Determinant = 0");
  S21Matrix buff;
  if (rows_ == 1 && cols_ == 1) {
    S21Matrix buff(*this);
  } else {
    S21Matrix buff = Transpose().CalcComplements();
  }
  buff.MulNumber(1 / det);
  return buff;
}

// overload

S21Matrix S21Matrix::operator+(const S21Matrix &A) const {
  S21Matrix sol(*this);
  sol.SumMatrix(A);
  return sol;
}

S21Matrix S21Matrix::operator-(const S21Matrix &A) const {
  S21Matrix sol(*this);
  sol.SubMatrix(A);
  return sol;
}

S21Matrix S21Matrix::operator*(const double number) const noexcept {
  S21Matrix sol(*this);
  sol.MulNumber(number);
  return sol;
}

S21Matrix S21Matrix::operator*(const S21Matrix &A) const {
  S21Matrix sol(*this);
  sol.MulMatrix(A);
  return sol;
}

bool S21Matrix::operator==(const S21Matrix &A) const noexcept {
  return EqMatrix(A);
}

bool S21Matrix::operator!=(const S21Matrix &A) const noexcept {
  return !(EqMatrix(A));
}

S21Matrix &S21Matrix::operator=(S21Matrix &A) {
  if (this != &A) {
    DeleteMatrix(*this);
    rows_ = A.rows_;
    cols_ = A.cols_;
    CopyMatrix(A);
  }
  return *this;
}

S21Matrix &S21Matrix::operator=(S21Matrix &&A) noexcept {
  if (this != &A) {
    DeleteMatrix(*this);
    rows_ = A.rows_;
    cols_ = A.cols_;
    matrix_ = A.matrix_;

    A.matrix_ = nullptr;
    A.rows_ = 0;
    A.cols_ = 0;
  }
  return *this;
}

S21Matrix S21Matrix::operator+=(const S21Matrix &A) {
  SumMatrix(A);
  return *this;
}

S21Matrix S21Matrix::operator-=(const S21Matrix &A) {
  SubMatrix(A);
  return *this;
}

S21Matrix S21Matrix::operator*=(const S21Matrix &A) {
  MulMatrix(A);
  return *this;
}

S21Matrix S21Matrix::operator*=(const double number) noexcept {
  MulNumber(number);
  return *this;
}

double &S21Matrix::operator()(int i, int j) {
  if (i >= rows_ || j >= cols_)
    throw std::out_of_range("Error! Value is out of range");
  if (i < 0 || j < 0)
    throw std::out_of_range("Error! Values should be positive");
  return matrix_[i][j];
}

const double &S21Matrix::operator()(int i, int j) const {
  if (i >= rows_ || j >= cols_)
    throw std::out_of_range("Error! Value is out of range");
  if (i < 0 || j < 0)
    throw std::out_of_range("Error! Values should be positive");
  return matrix_[i][j];
}

S21Matrix operator*(const double num, const S21Matrix &A) noexcept {
  S21Matrix sol(A);
  sol.MulNumber(num);
  return sol;
}
