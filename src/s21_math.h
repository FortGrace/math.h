#ifndef SRC_S21_MATH_H_
#define SRC_S21_MATH_H_

#define s21_PI 3.14159265358979323846264338327950288
#define s21_EXP 2.7182818284590452353602874713526624
#define s21_EPS 1e-99l
#define s21_NAN 0.0 / 0.0
#define s21_INF 1.0 / 0.0
#define s21_MAX_double 1.7976931348623158e308

#define is_fin(x) __builtin_isfinite(x)
#define is_nan(x) __builtin_isnan(x)
#define is_inf(x) __builtin_isinf(x)

int s21_abs(int x);
long double s21_ceil(double x);
long double s21_exp(double x);
long double s21_fabs(long double x);
long double s21_floor(double x);
long double s21_log(double x);
long double s21_pow(double base, double exp);
long double s21_sin(double x);
long double s21_cos(double x);
long double s21_sqrt(double x);
long double s21_tan(double x);
long double s21_asin(double x);
long double s21_acos(double x);
long double s21_atan(double x);
long double s21_fmod(double x, double y);

#endif  // SRC_S21_MATH_H_
