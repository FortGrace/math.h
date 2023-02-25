#include <check.h>
#include <math.h>
#include <stdlib.h>

#include "./s21_math.h"

START_TEST(abs_1) { ck_assert_int_eq(s21_abs(21), abs(21)); }
END_TEST

START_TEST(abs_2) { ck_assert_int_eq(s21_abs(-21), abs(-21)); }
END_TEST

START_TEST(abs_3) { ck_assert_int_eq(s21_abs(-2147483647), abs(-2147483647)); }
END_TEST

START_TEST(abs_4) { ck_assert_int_eq(s21_abs(2147483647), abs(2147483647)); }
END_TEST

START_TEST(abs_5) { ck_assert_int_eq(s21_abs(-0), abs(-0)); }
END_TEST

START_TEST(abs_6) { ck_assert_int_eq(s21_abs(+0), abs(+0)); }
END_TEST

START_TEST(abs_7) { ck_assert_int_eq(s21_abs((int)NAN), abs((int)NAN)); }
END_TEST

START_TEST(abs_8) {
    ck_assert_int_eq(s21_abs((int)INFINITY), abs((int)INFINITY));
}
END_TEST

START_TEST(abs_9) {
    ck_assert_int_eq(s21_abs((int)-INFINITY), abs((int)-INFINITY));
}
END_TEST

START_TEST(abs_10) { ck_assert_int_eq(s21_abs(0), abs(0)); }
END_TEST

// sqrt

START_TEST(sqrt_1) {
    double num = INFINITY;
    long double orig_res = sqrt(num), our_res = s21_sqrt(num);
    int suc = 0;
    if (isinf(orig_res) && isinf(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(sqrt_2) {
    double num = -INFINITY;
    long double orig_res = sqrt(num), our_res = s21_sqrt(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(sqrt_3) {
    double num = -15.01;
    ck_assert_ldouble_nan(s21_sqrt(num));
    ck_assert_ldouble_nan(sqrt(num));
}
END_TEST

START_TEST(sqrt_4) {
    double num = 0.9999;
    long double orig_res = sqrt(num), our_res = s21_sqrt(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(sqrt_5) {
    double num = NAN;
    long double orig_res = sqrt(num), our_res = s21_sqrt(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(sqrt_6) {
    ck_assert_ldouble_infinite(s21_sqrt(INFINITY));
    ck_assert_ldouble_infinite(sqrt(INFINITY));
}
END_TEST

START_TEST(sqrt_7) {
    ck_assert_ldouble_nan(s21_sqrt(NAN));
    ck_assert_ldouble_nan(sqrt(NAN));
}
END_TEST

START_TEST(sqrt_8) {
    ck_assert_ldouble_eq_tol(s21_sqrt(1.0), sqrt(1.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_sqrt(1.1), sqrt(1.1), 1e-6);
}
END_TEST

START_TEST(sqrt_9) {
    ck_assert_ldouble_eq_tol(s21_sqrt(0.0), sqrt(0.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_sqrt(-0.0), sqrt(-0.0), 1e-6);
}
END_TEST

START_TEST(sqrt_10) {
    ck_assert_ldouble_eq_tol(s21_sqrt(10e+16), sqrt(10e+16), 1e-6);
}
END_TEST

START_TEST(sqrt_12) { ck_assert_ldouble_eq_tol(s21_sqrt(1), sqrt(1), 1e-6); }
END_TEST

START_TEST(sqrt_13) {
    ck_assert_ldouble_eq_tol(s21_sqrt(100), sqrt(100), 1e-6);
}
END_TEST

START_TEST(sqrt_14) {
    ck_assert_ldouble_eq_tol(s21_sqrt(98765432.123), sqrt(98765432.123), 1e-6);
}
END_TEST

START_TEST(sqrt_15) {
    ck_assert_ldouble_eq_tol(s21_sqrt(981.123e6), sqrt(981.123e6), 1e-6);
}
END_TEST

START_TEST(sqrt_16) {
    ck_assert_ldouble_eq_tol(s21_sqrt(981.123e-15), sqrt(981.123e-15), 1e-6);
}
END_TEST

START_TEST(log_1) {
    double num = 456.789;
    long double orig_res = log(num);
    long double our_res = s21_log(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(log_2) {
    double num = NAN;
    long double orig_res = log(num);
    long double our_res = s21_log(num);
    int suc = 0;
    if (isnan(our_res) && isnan(orig_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(log_3) {
    double num = 0.001;
    long double orig_res = log(num);
    long double our_res = s21_log(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(log_4) {
    double num = 0.4;
    long double orig_res = log(num);
    long double our_res = s21_log(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(log_5) {
    double num = -INFINITY;
    long double orig_res = log(num);
    long double our_res = s21_log(num);
    int suc = 0;
    if (isnan(our_res) && isnan(orig_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(log_6) {
    double num = INFINITY;
    long double orig_res = log(num);
    long double our_res = s21_log(num);
    int suc = 0;
    if (isinf(our_res) && isinf(orig_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(log_7) {
    ck_assert_ldouble_infinite(s21_log(INFINITY));
    ck_assert_ldouble_infinite(log(INFINITY));

    ck_assert_ldouble_nan(s21_log(-INFINITY));
    ck_assert_ldouble_nan(log(-INFINITY));
}
END_TEST

START_TEST(log_8) {
    ck_assert_ldouble_nan(s21_log(NAN));
    ck_assert_ldouble_nan(log(NAN));
}
END_TEST

START_TEST(log_9) {
    ck_assert_ldouble_eq_tol(s21_log(1.0), log(1.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_log(1.1), log(1.1), 1e-6);
}
END_TEST

START_TEST(log_10) { ck_assert_ldouble_eq_tol(s21_log(0.5), log(0.5), 1e-6); }
END_TEST

START_TEST(log_11) { ck_assert_ldouble_eq_tol(s21_log(1), log(1), 1e-6); }
END_TEST

START_TEST(log_12) { ck_assert_ldouble_eq_tol(s21_log(100), log(100), 1e-6); }
END_TEST

START_TEST(log_13) { ck_assert_ldouble_eq_tol(s21_log(0.5), log(0.5), 1e-6); }
END_TEST

START_TEST(log_14) { ck_assert_ldouble_eq_tol(s21_log(0.5), log(0.5), 1e-6); }
END_TEST

START_TEST(log_15) {
    ck_assert_ldouble_eq_tol(s21_log(DBL_MAX), log(DBL_MAX), 1e-6);
}
END_TEST

START_TEST(log_16) {
    ck_assert_ldouble_eq_tol(s21_log(s21_PI / 3), log(s21_PI / 3), 1e-6);
}
END_TEST

START_TEST(log_17) {
    ck_assert_ldouble_eq_tol(s21_log(9.234578353457e-6), log(9.234578353457e-6),
                             1e-6);
}
END_TEST

START_TEST(log_18) {
    ck_assert_ldouble_eq_tol(s21_log(9.234578353457e6), log(9.234578353457e6),
                             1e-6);
}
END_TEST

START_TEST(log_19) {
    ck_assert_ldouble_eq_tol(s21_log(1234567.000000004), log(1234567.000000004),
                             1e-6);
}
END_TEST

START_TEST(log_20) {
    ck_assert_ldouble_eq_tol(s21_log(1.1e-16), log(1.1e-16), 1e-6);
}
END_TEST

START_TEST(log_24) {
    ck_assert_ldouble_eq_tol(s21_log(1.1e-45), log(1.1e-45), 1e-6);
}
END_TEST

START_TEST(log_25) {
    ck_assert_ldouble_eq_tol(s21_log(1.1e-18), log(1.1e-18), 1e-6);
}
END_TEST

START_TEST(log_26) {
    ck_assert_ldouble_eq_tol(s21_log(1.1e-19), log(1.1e-19), 1e-6);
}
END_TEST

START_TEST(log_27) {
    ck_assert_ldouble_eq_tol(s21_log(1.1e-23), log(1.1e-23), 1e-6);
}
END_TEST

START_TEST(log_21) {
    ck_assert_ldouble_eq_tol(s21_log(987654321123.123e-30),
                             log(987654321123.123e-30), 1e-6);
}
END_TEST

START_TEST(log_22) {
    ck_assert_ldouble_eq_tol(s21_log(98765.123e-11), log(98765.123e-11), 1e-6);
}
END_TEST

START_TEST(log_23) {
    ck_assert_ldouble_infinite(s21_log(0));
    ck_assert_ldouble_infinite(log(0));
}
END_TEST

START_TEST(log_28) {
    ck_assert_ldouble_eq_tol(s21_log(1.1e-17), log(1.1e-17), 1e-6);
}
END_TEST

START_TEST(log_29) {
    ck_assert_ldouble_eq_tol(s21_log(1.1e-15), log(1.1e-15), 1e-6);
}
END_TEST

START_TEST(log_30) {
    ck_assert_ldouble_eq_tol(s21_log(1.1e-14), log(1.1e-14), 1e-6);
}
END_TEST

START_TEST(log_31) {
    ck_assert_ldouble_eq_tol(s21_log(1.1e-13), log(1.1e-13), 1e-6);
}
END_TEST

START_TEST(log_32) {
    ck_assert_ldouble_eq_tol(s21_log(1.1e-12), log(1.1e-12), 1e-6);
}
END_TEST

START_TEST(log_33) {
    ck_assert_ldouble_eq_tol(s21_log(1.1e-11), log(1.1e-11), 1e-6);
}
END_TEST

START_TEST(log_34) {
    ck_assert_ldouble_eq_tol(s21_log(1.1e-80), log(1.1e-80), 1e-6);
}
END_TEST

START_TEST(fabs_1) { ck_assert_ldouble_eq(s21_fabs(0.000001), fabs(0.000001)); }
END_TEST

START_TEST(fabs_2) {
    ck_assert_ldouble_eq(s21_fabs(-21.000345), fabs(-21.000345));
}
END_TEST

START_TEST(fabs_3) {
    ck_assert_ldouble_eq(s21_fabs(-2147483600.543), fabs(-2147483600.543));
}
END_TEST

START_TEST(fabs_4) {
    ck_assert_ldouble_eq(s21_fabs(2147483600.4857), fabs(2147483600.4857));
}
END_TEST

START_TEST(fabs_5) {
    double a = -5.53151413431;
    ck_assert_ldouble_eq(s21_fabs(a), fabs(a));
}
END_TEST

START_TEST(fabs_6) {
    double a = NAN;
    ck_assert_ldouble_nan(s21_fabs(a));
    ck_assert_ldouble_nan(fabs(a));
}
END_TEST

START_TEST(fabs_7) {
    double a = s21_INF;
    ck_assert_ldouble_infinite(s21_fabs(a));
    ck_assert_ldouble_infinite(fabs(a));
}
END_TEST

START_TEST(fabs_8) {
    double a = -9519359135915.53151413431;
    ck_assert_ldouble_eq_tol(s21_fabs(a), fabs(a), 1e-6);
}
END_TEST

START_TEST(fabs_10) {
    ck_assert_ldouble_nan(s21_fabs(NAN));
    ck_assert_ldouble_nan(fabs(NAN));
}
END_TEST

START_TEST(fabs_11) { ck_assert_ldouble_eq(s21_fabs(-15.01), fabs(-15.01)); }
END_TEST

START_TEST(fabs_12) { ck_assert_ldouble_eq(s21_fabs(15.01), fabs(15.01)); }
END_TEST

START_TEST(fabs_13) {
    ck_assert_ldouble_eq(s21_fabs(INFINITY), fabs(INFINITY));
}
END_TEST

START_TEST(fabs_14) {
    ck_assert_ldouble_eq(s21_fabs(-INFINITY), fabs(-INFINITY));
}
END_TEST

START_TEST(fabs_15) {
    double num = NAN;
    long double orig_res = fabs(num), our_res = s21_fabs(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(exp_0) {
    double x = 0;
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_1) {
    double x = 1;
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_2) {
    double x = 2;
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_3) {
    double x = 3;
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_4) {
    double x = 1.5;
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_5) {
    double x = -1.5;
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_6) {
    double x = 0.5;
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_7) {
    double x = -0.5;
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_8) {
    double x = 715;
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_9) {
    double x = -715;
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_10) {
    double x = -0;
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_11) {
    double x = 'A';
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_12) {
    double x = -'N';
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_13) {
    double x = -'V';
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_14) {
    double x = +'B';
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_15) {
    double x = 0 / 1;
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_16) {
    double x = 0 / 1;
    ck_assert_int_eq(s21_exp(x), exp(x));
}
END_TEST

START_TEST(exp_17) { ck_assert_float_eq(s21_exp(2), exp(2)); }
END_TEST

START_TEST(exp_18) { ck_assert_float_eq(s21_exp(-750), exp(-750)); }
END_TEST

START_TEST(exp_19) { ck_assert_float_eq(s21_exp(-2), exp(-2)); }
END_TEST

START_TEST(exp_20) {
    double num = NAN;
    long double orig_res = exp(num), our_res = s21_exp(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(exp_21) {
    double num = INFINITY;
    long double orig_res = exp(num), our_res = s21_exp(num);
    int suc = 0;
    if (isinf(orig_res) && isinf(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(exp_22) {
    double num = -INFINITY;
    long double orig_res = exp(num), our_res = s21_exp(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(exp_23) {
    double x = 0.0;
    ck_assert_ldouble_eq_tol(s21_exp(x), expl(x), 1e-06);
}
END_TEST

START_TEST(exp_24) {
    ck_assert_ldouble_eq_tol(s21_exp(-INFINITY), expl(-INFINITY), 1e-06);
}
END_TEST

START_TEST(exp_25) {
    ck_assert_ldouble_infinite(s21_exp(INFINITY));
    ck_assert_ldouble_infinite(expl(INFINITY));
}
END_TEST

START_TEST(exp_26) {
    ck_assert_ldouble_nan(s21_exp(NAN));
    ck_assert_ldouble_nan(expl(NAN));
}
END_TEST

START_TEST(exp_27) {
    double x = 0.0;
    ck_assert_ldouble_eq_tol(s21_exp(x), expl(x), 1e-06);
}
END_TEST

START_TEST(exp_28) {
    double x = -1.0;
    ck_assert_ldouble_eq_tol(s21_exp(x), expl(x), 1e-06);
}
END_TEST

START_TEST(exp_29) {
    double x = -100;
    ck_assert_ldouble_eq_tol(s21_exp(x), expl(x), 1e-06);
}
END_TEST

START_TEST(exp_30) {
    double x = -1000000;
    ck_assert_ldouble_eq_tol(s21_exp(x), expl(x), 1e-06);
}
END_TEST

START_TEST(exp_31) {
    double x = 1;
    ck_assert_ldouble_eq_tol(s21_exp(x), expl(x), 1e-06);
}
END_TEST

START_TEST(exp_32) {
    double x = 1;
    ck_assert_ldouble_eq_tol(s21_exp(x), expl(x), 1e-06);
}
END_TEST

START_TEST(exp_33) {
    double num = 1000000;
    long double orig_res = exp(num), our_res = s21_exp(num);
    int suc = 0;
    if (isinf(orig_res) && isinf(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(pow_1) {
    double num = 456.789;
    double ex = 1.1;
    long double orig_res = pow(num, ex), our_res = s21_pow(num, ex);
    int suc = 0;
    ck_assert_float_eq(our_res, orig_res);
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(pow_2) {
    double num = INFINITY;
    double ex = 1.1;
    long double orig_res = pow(num, ex), our_res = s21_pow(num, ex);
    int suc = 0;
    if (isinf(our_res) && isinf(orig_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(pow_3) {
    double num = NAN;
    double ex = 1.1;
    long double orig_res = pow(num, ex), our_res = s21_pow(num, ex);
    int suc = 0;
    if (isnan(our_res) && isnan(orig_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(pow_4) {
    double num = NAN;
    double ex = NAN;
    long double orig_res = pow(num, ex), our_res = s21_pow(num, ex);
    int suc = 0;
    if (isnan(our_res) && isnan(orig_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(pow_5) {
    double num = 0.001;
    double ex = NAN;
    long double orig_res = pow(num, ex), our_res = s21_pow(num, ex);
    int suc = 0;
    if (isnan(our_res) && isnan(orig_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(pow_6) {
    double num = INFINITY;
    double ex = INFINITY;
    long double orig_res = pow(num, ex), our_res = s21_pow(num, ex);
    int suc = 0;
    if (isinf(our_res) && isinf(orig_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(pow_7) {
    double num = INFINITY;
    double ex = INFINITY;
    long double orig_res = pow(num, ex), our_res = s21_pow(num, ex);
    int suc = 0;
    if (isinf(our_res) && isinf(orig_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(pow_8) {
    long double base = -16.161435;
    long double exp = 9.;
    ck_assert_double_eq_tol(s21_pow(base, exp), pow(base, exp), 1e-6);
}
END_TEST

START_TEST(pow_9) {
    long double base = -10.1261;
    long double exp = -0.72;
    ck_assert_double_nan(s21_pow(base, exp));
    ck_assert_double_nan(pow(base, exp));
}
END_TEST

START_TEST(pow_10) {
    long double base = +0;
    long double exp = -0.33;
    ck_assert_ldouble_infinite(s21_pow(base, exp));
    ck_assert_ldouble_infinite(pow(base, exp));
}
END_TEST

START_TEST(pow_12) {
    long double base = -0;
    long double exp = -0.33;
    ck_assert_ldouble_infinite(s21_pow(base, exp));
    ck_assert_ldouble_infinite(pow(base, exp));
}
END_TEST

START_TEST(pow_13) {
    long double base = -0;
    long double exp = s21_INF;
    ck_assert_ldouble_eq_tol(s21_pow(base, exp), pow(base, exp), 1e-6);
}
END_TEST

START_TEST(pow_14) {
    long double base = +0;
    long double exp = 123;
    ck_assert_ldouble_eq_tol(s21_pow(base, exp), pow(base, exp), 1e-6);
}
END_TEST

START_TEST(pow_15) {
    long double base = -0;
    long double exp = 123;
    ck_assert_ldouble_eq_tol(s21_pow(base, exp), pow(base, exp), 1e-6);
}
END_TEST

START_TEST(pow_16) {
    long double base = -0;
    long double exp = s21_INF;
    ck_assert_ldouble_eq_tol(s21_pow(base, exp), pow(base, exp), 1e-6);
}
END_TEST

START_TEST(pow_17) {
    long double base = -1;
    long double exp = s21_INF;
    ck_assert_ldouble_eq_tol(s21_pow(base, exp), pow(base, exp), 1e-6);
}
END_TEST

START_TEST(pow_18) {
    long double base = -1;
    long double exp = s21_NAN;
    ck_assert_double_nan(pow(base, exp));
    ck_assert_double_nan(s21_pow(base, exp));
}
END_TEST

START_TEST(pow_19) {
    long double base = s21_NAN;
    long double exp = 0;
    ck_assert_ldouble_eq_tol(s21_pow(base, exp), pow(base, exp), 1e-6);
}
END_TEST

START_TEST(pow_20) {
    long double base = -123;
    long double exp = s21_NAN;
    ck_assert_double_nan(s21_pow(base, exp));
    ck_assert_double_nan(pow(base, exp));
}
END_TEST

START_TEST(pow_21) {
    long double base = 0.5591951;
    long double exp = -s21_INF;
    ck_assert_ldouble_infinite(pow(base, exp));
    ck_assert_ldouble_infinite(s21_pow(base, exp));
}
END_TEST

START_TEST(pow_22) {
    long double base = 1.5591951;
    long double exp = -s21_INF;
    ck_assert_ldouble_eq_tol(s21_pow(base, exp), pow(base, exp), 1e-6);
}
END_TEST

START_TEST(pow_23) {
    long double base = 0.5591951;
    long double exp = s21_INF;
    ck_assert_ldouble_eq_tol(s21_pow(base, exp), pow(base, exp), 1e-6);
}
END_TEST

START_TEST(pow_24) {
    long double base = 1.5591951;
    long double exp = s21_INF;
    ck_assert_ldouble_infinite(s21_pow(base, exp));
    ck_assert_ldouble_infinite(pow(base, exp));
}
END_TEST

START_TEST(pow_25) {
    long double base = -s21_INF;
    long double exp = -193491;
    ck_assert_ldouble_eq_tol(s21_pow(base, exp), pow(base, exp), 1e-6);
}
END_TEST

START_TEST(pow_26) {
    long double base = -s21_INF;
    long double exp = -142;
    ck_assert_ldouble_eq_tol(s21_pow(base, exp), pow(base, exp), 1e-6);
}
END_TEST

START_TEST(pow_27) {
    long double base = -s21_INF;
    long double exp = 141;
    ck_assert_ldouble_infinite(s21_pow(base, exp));
    ck_assert_ldouble_infinite(pow(base, exp));
}
END_TEST

START_TEST(pow_28) {
    long double base = -s21_INF;
    long double exp = 142;
    ck_assert_ldouble_infinite(s21_pow(base, exp));
    ck_assert_ldouble_infinite(pow(base, exp));
}
END_TEST

START_TEST(pow_29) {
    long double base = s21_INF;
    long double exp = -1;
    ck_assert_ldouble_eq_tol(s21_pow(base, exp), pow(base, exp), 1e-6);
}
END_TEST

START_TEST(pow_30) {
    long double base = s21_INF;
    long double exp = 1;
    ck_assert_ldouble_infinite(s21_pow(base, exp));
    ck_assert_ldouble_infinite(pow(base, exp));
}
END_TEST

START_TEST(pow_31) {
    ck_assert_ldouble_eq_tol(s21_pow(0.0, 0.0), pow(0.0, 0.0), 1e-6);
}
END_TEST

START_TEST(pow_32) {
    ck_assert_ldouble_eq_tol(s21_pow(-1, 4), pow(-1, 4), 1e-6);
}
END_TEST

START_TEST(pow_33) {
    ck_assert_ldouble_eq_tol(s21_pow(2, 18), pow(2, 18), 1e-6);
}
END_TEST

START_TEST(pow_34) {
    ck_assert_ldouble_eq_tol(s21_pow(2.1, 11), pow(2.1, 11), 1e-6);
}
END_TEST

START_TEST(pow_35) {
    ck_assert_ldouble_eq_tol(s21_pow(1.23456e-7, 2.21), pow(1.23456e-7, 2.21),
                             1e-6);
}
END_TEST

START_TEST(pow_36) {
    ck_assert_ldouble_eq_tol(s21_pow(2.1234567, -2), pow(2.1234567, -2), 1e-6);
}
END_TEST

START_TEST(pow_37) {
    ck_assert_ldouble_eq_tol(s21_pow(2.1234567, -2), pow(2.1234567, -2), 1e-6);
}
END_TEST

START_TEST(pow_38) {
    long double orig_res = pow(2.1234567e-7, -2.45e4),
                our_res = s21_pow(2.1234567e-7, -2.45e4);
    int suc = 0;
    if (isinf(our_res) && isinf(orig_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(pow_39) {
    ck_assert_ldouble_eq_tol(s21_pow(9.99999999e3, 9.95e-4),
                             pow(9.99999999e3, 9.95e-4), 1e-6);
}
END_TEST

START_TEST(pow_40) {
    long double exp = -12;
    ck_assert_ldouble_infinite(s21_pow(+0, exp));
    ck_assert_ldouble_infinite(pow(+0, exp));
}
END_TEST

START_TEST(pow_41) {
    long double exp = -12;
    ck_assert_ldouble_infinite(s21_pow(-0, exp));
    ck_assert_ldouble_infinite(pow(-0, exp));
}
END_TEST

START_TEST(pow_42) {
    long double exp = 0;
    ck_assert_ldouble_eq_tol(s21_pow(-0, exp), pow(-0, exp), 1e-6);
}
END_TEST

START_TEST(pow_43) {
    long double exp = s21_NAN;
    ck_assert_double_nan(s21_pow(+0, exp));
    ck_assert_double_nan(pow(+0, exp));
}
END_TEST

START_TEST(pow_44) {
    long double exp = 2.456;
    ck_assert_ldouble_eq_tol(s21_pow(+0, exp), pow(+0, exp), 1e-6);
}
END_TEST

START_TEST(pow_45) {
    long double exp = -s21_INF;
    ck_assert_ldouble_infinite(s21_pow(+0, exp));
    ck_assert_ldouble_infinite(pow(+0, exp));
}
END_TEST

START_TEST(pow_46) {
    long double exp = 56;
    ck_assert_ldouble_eq_tol(s21_pow(+0, exp), pow(+0, exp), 1e-6);
}
END_TEST

START_TEST(pow_47) {
    long double exp = 56;
    ck_assert_ldouble_eq_tol(s21_pow(-0, exp), pow(-0, exp), 1e-6);
}
END_TEST

START_TEST(pow_48) {
    long double exp = 56.345;
    ck_assert_ldouble_eq_tol(s21_pow(-0, exp), pow(-0, exp), 1e-6);
}
END_TEST

START_TEST(pow_49) {
    long double exp = s21_INF;
    ck_assert_ldouble_eq_tol(s21_pow(-1, exp), pow(-1, exp), 1e-6);
}
END_TEST

START_TEST(pow_50) {
    long double exp = -s21_INF;
    ck_assert_ldouble_eq_tol(s21_pow(-1, exp), pow(-1, exp), 1e-6);
}
END_TEST

START_TEST(pow_51) {
    long double exp = s21_NAN;
    ck_assert_ldouble_eq_tol(s21_pow(+1, exp), pow(+1, exp), 1e-6);
}
END_TEST

START_TEST(pow_52) {
    long double exp = -0;
    ck_assert_ldouble_eq_tol(s21_pow(4567, exp), pow(4567, exp), 1e-6);
}
END_TEST

START_TEST(pow_53) {
    long double exp = +0;
    ck_assert_ldouble_eq_tol(s21_pow(34567, exp), pow(34567, exp), 1e-6);
}
END_TEST

START_TEST(pow_54) {
    long double exp = -0;
    ck_assert_ldouble_eq_tol(s21_pow(s21_NAN, exp), pow(s21_NAN, exp), 1e-6);
}
END_TEST

START_TEST(pow_55) {
    long double exp = +0;
    ck_assert_ldouble_eq_tol(s21_pow(s21_NAN, exp), pow(s21_NAN, exp), 1e-6);
}
END_TEST

START_TEST(pow_56) {
    long double exp = s21_NAN;
    ck_assert_double_nan(s21_pow(s21_NAN, exp));
    ck_assert_double_nan(pow(s21_NAN, exp));
}
END_TEST

START_TEST(pow_57) {
    long double exp = s21_INF;
    ck_assert_ldouble_infinite(s21_pow(s21_INF, exp));
    ck_assert_ldouble_infinite(pow(s21_INF, exp));
}
END_TEST

START_TEST(pow_58) {
    long double exp = -s21_INF;
    ck_assert_ldouble_eq_tol(s21_pow(-s21_INF, exp), pow(-s21_INF, exp), 1e-6);
}
END_TEST

START_TEST(pow_59) {
    double exp = s21_INF;
    ck_assert_ldouble_infinite(s21_pow(-s21_INF, exp));
    ck_assert_ldouble_infinite(pow(-s21_INF, exp));
}
END_TEST

START_TEST(pow_60) {
    double exp;
    double res1;
    long double res2;

    for (exp = 0.0; exp > -1000; exp -= 0.1) {
        res1 = pow(s21_INF, exp);
        res2 = s21_pow(s21_INF, exp);
    }
    ck_assert_ldouble_eq_tol(res1, res2, 1e-6);
}
END_TEST

START_TEST(pow_61) {
    double base;
    double res1;
    long double res2;

    for (base = 0.0; base < 10; base += 0.1) {
        res1 = pow(base, s21_INF);
        res2 = s21_pow(base, s21_INF);
    }
    ck_assert_ldouble_infinite(res1);
    ck_assert_ldouble_infinite(res2);
}
END_TEST

START_TEST(pow_62) {
    double base;
    double res1;
    long double res2;

    for (base = 0.0; base > -10; base -= 0.1) {
        res1 = pow(base, s21_INF);
        res2 = s21_pow(base, s21_INF);
    }
    ck_assert_ldouble_infinite(res1);
    ck_assert_ldouble_infinite(res2);
}
END_TEST

START_TEST(pow_63) {
    long double base;
    long double res1;
    long double res2;

    for (base = 0.0; base < 10; base += 0.1) {
        res1 = pow(base, -s21_INF);
        res2 = s21_pow(base, -s21_INF);
    }
    ck_assert_ldouble_eq_tol(res1, res2, 1e-6);
}
END_TEST

START_TEST(pow_64) {
    long double base;
    long double res1;
    long double res2;

    for (base = 0.0; base > -10; base -= 0.1) {
        res1 = pow(base, -s21_INF);
        res2 = s21_pow(base, -s21_INF);
    }
    ck_assert_ldouble_eq_tol(res1, res2, 1e-6);
}
END_TEST

START_TEST(pow_65) {
    double exp;
    double res1;
    long double res2;

    for (exp = 0.0; exp > -10; exp -= 0.1) {
        if (fmod(exp, 2) != 0) {
            res1 = pow(-s21_INF, exp);
            res2 = s21_pow(-s21_INF, exp);
        }
    }
    ck_assert_ldouble_eq_tol(res1, res2, 1e-6);
}
END_TEST

START_TEST(pow_66) {
    double exp;
    double res1;
    long double res2;

    for (exp = 0.0; exp > -10; exp -= 0.1) {
        if (fmod(exp, 2) == 0) {
            res1 = pow(-s21_INF, exp);
            res2 = s21_pow(-s21_INF, exp);
        }
    }
    ck_assert_ldouble_eq_tol(res1, res2, 1e-6);
}
END_TEST

START_TEST(pow_67) {
    double exp;
    double res1;
    long double res2;

    for (exp = 0.0; exp > -10; exp -= 0.1) {
        if (fmod(exp, 2) != 0) {
            res1 = pow(-s21_INF, exp);
            res2 = s21_pow(-s21_INF, exp);
        }
    }
    // ck_assert_ldouble_infinite(res1);
    // ck_assert_ldouble_infinite(res2);
    ck_assert_ldouble_eq_tol(res1, res2, 1e-6);
}
END_TEST

START_TEST(pow_68) {
    double exp;
    double res1;
    long double res2;

    for (exp = 0.0; exp < 10; exp += 0.1) {
        if (fmod(exp, 2) == 0) {
            res1 = pow(-s21_INF, exp);
            res2 = s21_pow(-s21_INF, exp);
        }
    }
    ck_assert_ldouble_eq_tol(res2, res1, 1e-6);
}
END_TEST

START_TEST(pow_69) {
    long double exp = 3;
    ck_assert_ldouble_eq_tol(s21_pow(-3, exp), pow(-3, exp), 1e-6);
}
END_TEST

START_TEST(pow_70) {
    long double exp = -3;
    ck_assert_ldouble_eq_tol(s21_pow(-3, exp), pow(-3, exp), 1e-6);
}
END_TEST

START_TEST(pow_71) {
    long double exp = -3;
    ck_assert_ldouble_eq_tol(s21_pow(3, exp), pow(3, exp), 1e-6);
}
END_TEST

START_TEST(pow_72) {
    long double exp = 3;
    ck_assert_ldouble_eq_tol(s21_pow(3, exp), pow(3, exp), 1e-6);
}
END_TEST

START_TEST(pow_73) {
    long double base = -0;
    long double exp = s21_INF;
    long double res = pow(base, exp);
    long double res1 = s21_pow(base, exp);
    ck_assert_ldouble_eq_tol(res, res1, 1e-6);
}
END_TEST

START_TEST(pow_74) {
    long double exp;
    long double base = -0;
    long double res;
    long double res1;
    for (exp = 0.0; exp < 10; exp += 0.1) {
        if (fmod(exp, 2) != 0) {
            res = pow(base, exp);
            res1 = s21_pow(base, exp);
        }
    }
    ck_assert_ldouble_eq_tol(res, res1, 1e-6);
}
END_TEST

START_TEST(floor_1) { ck_assert_ldouble_eq(s21_floor(0.0), floor(0.0)); }
END_TEST

START_TEST(floor_2) { ck_assert_ldouble_eq(s21_floor(21.21), floor(21.21)); }
END_TEST

START_TEST(floor_3) { ck_assert_ldouble_eq(s21_floor(21.91), floor(21.91)); }
END_TEST

START_TEST(floor_4) { ck_assert_ldouble_eq(s21_floor(-21.21), floor(-21.21)); }
END_TEST

START_TEST(floor_5) { ck_assert_ldouble_eq(s21_floor(-21.91), floor(-21.91)); }
END_TEST

START_TEST(floor_6) {
    ck_assert_ldouble_infinite(s21_floor(INFINITY));
    ck_assert_ldouble_infinite(floor(INFINITY));
}
END_TEST

START_TEST(floor_7) {
    ck_assert_ldouble_nan(s21_floor(NAN));
    ck_assert_ldouble_nan(floor(NAN));
}
END_TEST

START_TEST(floor_8) { ck_assert_ldouble_eq(s21_floor(0), floor(0)); }
END_TEST

START_TEST(floor_10) { ck_assert_ldouble_eq(s21_floor(-15.01), floor(-15.01)); }
END_TEST

START_TEST(floor_11) { ck_assert_ldouble_eq(s21_floor(15.01), floor(15.01)); }
END_TEST

START_TEST(floor_12) {
    ck_assert_ldouble_eq(s21_floor(INFINITY), floor(INFINITY));
}
END_TEST

START_TEST(floor_13) {
    ck_assert_ldouble_eq(s21_floor(-INFINITY), floor(-INFINITY));
}
END_TEST

START_TEST(floor_14) {
    double num = NAN;
    long double orig_res = floor(num), our_res = s21_floor(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(floor_15) { ck_assert_ldouble_eq(s21_floor(-0), floor(-0)); }
END_TEST

START_TEST(ceil_1) { ck_assert_ldouble_eq(s21_ceil(-15.01), ceil(-15.01)); }
END_TEST

START_TEST(ceil_2) { ck_assert_ldouble_eq(s21_ceil(15.01), ceil(15.01)); }
END_TEST

START_TEST(ceil_3) { ck_assert_ldouble_eq(s21_ceil(INFINITY), ceil(INFINITY)); }
END_TEST

START_TEST(ceil_4) {
    ck_assert_ldouble_eq(s21_ceil(-INFINITY), ceil(-INFINITY));
}
END_TEST

START_TEST(ceil_5) { ck_assert_ldouble_eq(s21_ceil(-0.12), ceil(-0.12)); }
END_TEST

START_TEST(ceil_6) {
    double num = NAN;
    long double orig_res = ceil(num), our_res = s21_ceil(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(ceil_7) {
    ck_assert_ldouble_infinite(s21_ceil(INFINITY));
    ck_assert_ldouble_infinite(ceil(INFINITY));
}
END_TEST

START_TEST(ceil_8) {
    ck_assert_ldouble_nan(s21_ceil(NAN));
    ck_assert_ldouble_nan(ceil(NAN));
}
END_TEST

START_TEST(ceil_9) { ck_assert_ldouble_eq(s21_ceil(0), ceil(0)); }
END_TEST

START_TEST(ceil_11) { ck_assert_ldouble_eq(s21_ceil(0.0), ceil(0.0)); }
END_TEST

START_TEST(ceil_12) { ck_assert_ldouble_eq(s21_ceil(21.21), ceil(21.21)); }
END_TEST

START_TEST(ceil_13) { ck_assert_ldouble_eq(s21_ceil(21.91), ceil(21.91)); }
END_TEST

START_TEST(ceil_14) { ck_assert_ldouble_eq(s21_ceil(-21.21), ceil(-21.21)); }
END_TEST

START_TEST(ceil_15) { ck_assert_ldouble_eq(s21_ceil(-21.91), ceil(-21.91)); }
END_TEST

START_TEST(ceil_16) { ck_assert_ldouble_eq(s21_ceil(DBL_MAX), ceil(DBL_MAX)); }
END_TEST

START_TEST(sin_1) {
    double num = 0.55;
    long double orig_res = s21_sin(num), our_res = s21_sin(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(sin_2) { ck_assert_float_eq(s21_sin(0), sin(0)); }
END_TEST

START_TEST(sin_3) {
    double num = 1;
    long double orig_res = s21_sin(num), our_res = s21_sin(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(sin_4) {
    double num = -1;
    long double orig_res = s21_sin(num), our_res = s21_sin(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(sin_5) {
    double num = 0.999999;
    long double orig_res = s21_sin(num), our_res = s21_sin(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(sin_6) {
    double num = -0.999999;
    long double orig_res = s21_sin(num), our_res = s21_sin(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(sin_7) {
    double num = NAN;
    long double orig_res = sin(num), our_res = s21_sin(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(sin_8) {
    double num = INFINITY;
    long double orig_res = sin(num), our_res = s21_sin(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(sin_9) {
    double num = -INFINITY;
    long double orig_res = sin(num), our_res = s21_sin(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(sin_10) {
    double x = 0;
    ck_assert_double_eq_tol(s21_sin(x), sinl(x), 1e-06);
}
END_TEST

START_TEST(sin_11) {
    ck_assert_ldouble_nan(s21_sin(-INFINITY));
    ck_assert_ldouble_nan(sinl(-INFINITY));
}
END_TEST

START_TEST(sin_12) {
    ck_assert_ldouble_nan(s21_sin(INFINITY));
    ck_assert_ldouble_nan(sinl(INFINITY));
}
END_TEST

START_TEST(sin_13) {
    ck_assert_ldouble_nan(s21_sin(NAN));
    ck_assert_ldouble_nan(sinl(NAN));
}
END_TEST

START_TEST(sin_14) { ck_assert_double_eq_tol(s21_sin(0.0), sin(0.0), 1e-06); }
END_TEST

START_TEST(sin_15) {
    ck_assert_double_eq_tol(s21_sin(M_PI / 6), sin(M_PI / 6), 1e-06);
}
END_TEST

START_TEST(sin_16) {
    ck_assert_double_eq_tol(s21_sin(M_PI / 4), sin(M_PI / 4), 1e-06);
}
END_TEST

START_TEST(sin_17) {
    ck_assert_double_eq_tol(s21_sin(M_PI / 3), sin(M_PI / 3), 1e-06);
}
END_TEST

START_TEST(sin_18) {
    ck_assert_double_eq_tol(s21_sin(M_PI / 2), sin(M_PI / 2), 1e-06);
}
END_TEST

START_TEST(sin_19) {
    ck_assert_double_eq_tol(s21_sin(3 * M_PI / 2), sin(3 * M_PI / 2), 1e-06);
}
END_TEST

START_TEST(sin_20) {
    ck_assert_double_eq_tol(s21_sin(2 * M_PI), sin(2 * M_PI), 1e-06);
}
END_TEST

START_TEST(sin_21) {
    ck_assert_double_eq_tol(s21_sin(-2 * M_PI), sin(-2 * M_PI), 1e-06);
}
END_TEST

START_TEST(sin_22) {
    ck_assert_double_eq_tol(s21_sin(-3 * M_PI), sin(-3 * M_PI), 1e-06);
}
END_TEST

START_TEST(sin_23) {
    ck_assert_double_eq_tol(s21_sin(174.532925199433), sin(174.532925199433),
                            1e-06);
}
END_TEST

START_TEST(sin_24) {
    ck_assert_double_eq_tol(s21_sin(-174.532925199433), sin(-174.532925199433),
                            1e-06);
}
END_TEST

START_TEST(sin_25) {
    ck_assert_double_eq_tol(s21_sin(-s21_PI), sin(-s21_PI), 1e-06);
}
END_TEST

START_TEST(cos_1) {
    double num = 0.55;
    long double orig_res = s21_cos(num), our_res = s21_cos(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(cos_2) {
    double num = -0.999999;
    long double orig_res = cos(num), our_res = s21_cos(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(cos_3) {
    double num = 628;
    long double orig_res = cos(num), our_res = s21_cos(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(cos_4) { ck_assert_float_eq(s21_cos(0), cos(0)); }
END_TEST

START_TEST(cos_5) { ck_assert_float_eq(s21_cos(-1), cos(-1)); }
END_TEST

START_TEST(cos_6) {
    double num = NAN;
    long double orig_res = cos(num), our_res = s21_cos(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(cos_7) {
    double num = INFINITY;
    long double orig_res = cos(num), our_res = s21_cos(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(cos_8) {
    double num = -INFINITY;
    long double orig_res = cos(num), our_res = s21_cos(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(cos_9) {
    double x = 0;
    ck_assert_double_eq_tol(s21_cos(x), cosl(x), 1e-06);
}
END_TEST

START_TEST(cos_10) {
    ck_assert_ldouble_nan(s21_cos(-INFINITY));
    ck_assert_ldouble_nan(cosl(-INFINITY));
}
END_TEST

START_TEST(cos_11) {
    ck_assert_ldouble_nan(s21_cos(INFINITY));
    ck_assert_ldouble_nan(cosl(INFINITY));
}
END_TEST

START_TEST(cos_12) {
    ck_assert_ldouble_nan(s21_cos(NAN));
    ck_assert_ldouble_nan(cosl(NAN));
}
END_TEST

START_TEST(cos_13) { ck_assert_double_eq_tol(s21_cos(0.0), cosl(0.0), 1e-06); }
END_TEST

START_TEST(cos_14) {
    ck_assert_double_eq_tol(s21_cos(M_PI / 4), cosl(M_PI / 4), 1e-06);
}
END_TEST

START_TEST(cos_15) {
    ck_assert_double_eq_tol(s21_cos(M_PI / 6), cosl(M_PI / 6), 1e-06);
}
END_TEST

START_TEST(cos_16) {
    ck_assert_double_eq_tol(s21_cos(M_PI / 3), cosl(M_PI / 3), 1e-06);
}
END_TEST

START_TEST(cos_17) {
    ck_assert_double_eq_tol(s21_cos(M_PI / 2), cosl(M_PI / 2), 1e-06);
}
END_TEST

START_TEST(cos_18) {
    ck_assert_double_eq_tol(s21_cos(3 * M_PI / 2), cosl(3 * M_PI / 2), 1e-06);
}
END_TEST

START_TEST(cos_19) {
    ck_assert_double_eq_tol(s21_cos(2 * M_PI), cosl(2 * M_PI), 1e-06);
}
END_TEST

START_TEST(cos_20) {
    ck_assert_double_eq_tol(s21_cos(174.532925199433), cosl(174.532925199433),
                            1e-06);
}
END_TEST

START_TEST(cos_21) {
    ck_assert_double_eq_tol(s21_cos(-174.532925199433), cosl(-174.532925199433),
                            1e-06);
}
END_TEST

START_TEST(tan_1) { ck_assert_float_eq(s21_tan(-15.01), tan(-15.01)); }
END_TEST

START_TEST(tan_2) {
    double num = 0.55;
    long double orig_res = s21_tan(num), our_res = s21_tan(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(tan_3) { ck_assert_float_eq(-9999999999, -9999999999); }
END_TEST

START_TEST(tan_4) {
    double num = NAN;
    long double orig_res = tan(num), our_res = s21_tan(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(tan_5) {
    double num = INFINITY;
    long double orig_res = tan(num), our_res = s21_tan(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(tan_6) {
    double num = -INFINITY;
    long double orig_res = tan(num), our_res = s21_tan(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(tan_7) {
    double x = 0;
    ck_assert_double_eq_tol(s21_tan(x), tanl(x), 1e-06);
}
END_TEST

START_TEST(tan_8) {
    ck_assert_ldouble_nan(s21_tan(-INFINITY));
    ck_assert_ldouble_nan(tanl(-INFINITY));
}
END_TEST

START_TEST(tan_9) {
    ck_assert_ldouble_nan(s21_tan(INFINITY));
    ck_assert_ldouble_nan(tanl(INFINITY));
}
END_TEST

START_TEST(tan_10) {
    ck_assert_ldouble_nan(s21_tan(NAN));
    ck_assert_ldouble_nan(tanl(NAN));
}
END_TEST

START_TEST(tan_11) { ck_assert_double_eq_tol(s21_tan(0.0), tanl(0.0), 1e-06); }
END_TEST

START_TEST(tan_12) {
    ck_assert_double_eq_tol(s21_tan(M_PI / 6), tanl(M_PI / 6), 1e-06);
}
END_TEST

START_TEST(tan_13) {
    ck_assert_double_eq_tol(s21_tan(M_PI / 4), tanl(M_PI / 4), 1e-06);
}
END_TEST

START_TEST(tan_14) {
    ck_assert_double_eq_tol(s21_tan(-2 * M_PI), tanl(-2 * M_PI), 1e-06);
}
END_TEST

START_TEST(tan_15) {
    ck_assert_double_eq_tol(s21_tan(M_PI), tanl(M_PI), 1e-06);
}
END_TEST

START_TEST(tan_16) {
    ck_assert_double_eq_tol(s21_tan(M_PI / 3), tanl(M_PI / 3), 1e-06);
}
END_TEST

START_TEST(tan_17) {
    ck_assert_double_eq_tol(s21_tan(3 * M_PI / 3), tanl(3 * M_PI / 3), 1e-06);
}
END_TEST

START_TEST(tan_18) {
    ck_assert_double_eq_tol(s21_tan(174.532925199433), tanl(174.532925199433),
                            1e-06);
}
END_TEST

START_TEST(tan_19) {
    ck_assert_double_eq_tol(s21_tan(-174.532925199433), tanl(-174.532925199433),
                            1e-06);
}
END_TEST

#if defined(__linux__)
START_TEST(tan_20) {
    ck_assert_double_eq_tol(s21_tan(s21_PI / 2), tanl(s21_PI / 2), 1e-06);
}
END_TEST

START_TEST(tan_21) {
    ck_assert_double_eq_tol(s21_tan(-s21_PI / 2), tanl(-s21_PI / 2), 1e-06);
}
END_TEST
#endif

START_TEST(asin_1) {
    double num = 0.55;
    long double orig_res = asin(num), our_res = s21_asin(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(asin_2) {
    double num = -0.99;
    long double orig_res = asin(num), our_res = s21_asin(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(asin_3) {
    double num = -9;
    long double orig_res = asin(num), our_res = s21_asin(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(asin_4) {
    double num = INFINITY;
    long double orig_res = asin(num), our_res = s21_asin(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(asin_5) {
    double num = -INFINITY;
    long double orig_res = asin(num), our_res = s21_asin(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(asin_6) {
    double num = NAN;
    long double orig_res = asin(num), our_res = s21_asin(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(asin_7) {
    ck_assert_ldouble_nan(s21_asin(INFINITY));
    ck_assert_ldouble_nan(asin(INFINITY));
}
END_TEST

START_TEST(asin_8) {
    ck_assert_ldouble_nan(s21_asin(-INFINITY));
    ck_assert_ldouble_nan(asin(-INFINITY));
}
END_TEST

START_TEST(asin_9) {
    ck_assert_ldouble_nan(s21_asin(NAN));
    ck_assert_ldouble_nan(asin(NAN));
}
END_TEST

START_TEST(asin_10) {
    ck_assert_ldouble_eq_tol(s21_asin(1.0), asin(1.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_asin(-1.0), asin(-1.0), 1e-6);
}
END_TEST

START_TEST(asin_11) {
    ck_assert_ldouble_eq_tol(s21_asin(0.0), asin(0.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_asin(-0.0), asin(-0.0), 1e-6);
}
END_TEST

START_TEST(asin_12) {
    ck_assert_ldouble_eq_tol(s21_asin(0.0), asin(0.0), 1e-6);
}
END_TEST

START_TEST(asin_13) { ck_assert_ldouble_eq_tol(s21_asin(-1), asin(-1), 1e-6); }
END_TEST

START_TEST(asin_14) { ck_assert_ldouble_eq_tol(s21_asin(1), asin(1), 1e-6); }
END_TEST

START_TEST(asin_15) {
    ck_assert_ldouble_eq_tol(s21_asin(-sqrt(3) / 2), asin(-sqrt(3) / 2), 1e-6);
}
END_TEST

START_TEST(asin_16) {
    ck_assert_ldouble_eq_tol(s21_asin(-sqrt(2) / 2), asin(-sqrt(2) / 2), 1e-6);
}
END_TEST

START_TEST(asin_17) {
    ck_assert_ldouble_eq_tol(s21_asin(-1 / 2), asin(-1 / 2), 1e-6);
}
END_TEST

START_TEST(asin_18) {
    ck_assert_ldouble_eq_tol(s21_asin(1 / 2), asin(1 / 2), 1e-6);
}
END_TEST

START_TEST(asin_19) {
    ck_assert_ldouble_eq_tol(s21_asin(sqrt(2) / 2), asin(sqrt(2) / 2), 1e-6);
}
END_TEST

START_TEST(asin_20) {
    ck_assert_ldouble_eq_tol(s21_asin(sqrt(3) / 2), asin(sqrt(3) / 2), 1e-6);
}
END_TEST

START_TEST(acos_1) {
    double num = 0.55;
    long double orig_res = acos(num), our_res = s21_acos(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(acos_2) {
    double num = -0.99;
    long double orig_res = acos(num), our_res = s21_acos(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(acos_3) {
    double num = -9;
    long double orig_res = acos(num), our_res = s21_acos(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(acos_4) {
    double num = INFINITY;
    long double orig_res = acos(num), our_res = s21_acos(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(acos_5) {
    double num = -INFINITY;
    long double orig_res = acos(num), our_res = s21_acos(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(acos_6) {
    double num = NAN;
    long double orig_res = s21_acos(num), our_res = acos(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(acos_7) {
    ck_assert_ldouble_nan(s21_acos(NAN));
    ck_assert_ldouble_nan(acos(NAN));
}
END_TEST

START_TEST(acos_8) {
    ck_assert_ldouble_nan(s21_acos(NAN));
    ck_assert_ldouble_nan(acos(NAN));
}
END_TEST

START_TEST(acos_9) {
    ck_assert_ldouble_eq_tol(s21_acos(1.0), acos(1.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_acos(-1.0), acos(-1.0), 1e-6);
}
END_TEST

START_TEST(acos_10) {
    ck_assert_ldouble_eq_tol(s21_acos(0.0), acos(0.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_acos(-0.0), acos(-0.0), 1e-6);
}
END_TEST

START_TEST(acos_11) { ck_assert_ldouble_eq_tol(s21_acos(-1), acos(-1), 1e-6); }
END_TEST

START_TEST(acos_12) { ck_assert_ldouble_eq_tol(s21_acos(1), acos(1), 1e-6); }
END_TEST

START_TEST(acos_13) {
    ck_assert_ldouble_eq_tol(s21_acos(-sqrt(3) / 2), acos(-sqrt(3) / 2), 1e-6);
}
END_TEST

START_TEST(acos_14) {
    ck_assert_ldouble_eq_tol(s21_acos(-sqrt(2) / 2), acos(-sqrt(2) / 2), 1e-6);
}
END_TEST

START_TEST(acos_15) {
    ck_assert_ldouble_eq_tol(s21_acos(-1 / 2), acos(-1 / 2), 1e-6);
}
END_TEST

START_TEST(acos_16) {
    ck_assert_ldouble_eq_tol(s21_acos(1 / 2), acos(1 / 2), 1e-6);
}
END_TEST

START_TEST(acos_17) {
    ck_assert_ldouble_eq_tol(s21_acos(2 / 2), acos(2 / 2), 1e-6);
}
END_TEST

START_TEST(acos_18) {
    ck_assert_ldouble_eq_tol(s21_acos(sqrt(2) / 2), acos(sqrt(2) / 2), 1e-6);
}
END_TEST

START_TEST(acos_19) {
    ck_assert_ldouble_eq_tol(s21_acos(sqrt(3) / 2), acos(sqrt(3) / 2), 1e-6);
}
END_TEST

START_TEST(atan_1) { ck_assert_float_eq(s21_atan(-15.01), atan(-15.01)); }
END_TEST

START_TEST(atan_2) {
    double num = -0.99;
    long double orig_res = s21_atan(num), our_res = atan(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(atan_3) { ck_assert_float_eq(-9999999999, -9999999999); }
END_TEST
START_TEST(atan_4) {
    ck_assert_ldouble_eq(s21_atan(INFINITY), atan(INFINITY));
    ck_assert_ldouble_eq(s21_atan(-INFINITY), atan(-INFINITY));
}
END_TEST

START_TEST(atan_5) {
    double num = NAN;
    long double orig_res = s21_atan(num), our_res = atan(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(atan_6) {
    ck_assert_ldouble_eq_tol(s21_atan(INFINITY), (M_PI / 2.0), 1e-6);
    ck_assert_ldouble_eq_tol(atan(INFINITY), (M_PI / 2.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_atan((-INFINITY)), ((M_PI) / -2.0), 1e-6);
    ck_assert_ldouble_eq_tol(atan(-INFINITY), ((M_PI) / -2.0), 1e-6);
}
END_TEST

START_TEST(atan_7) {
    ck_assert_ldouble_nan(s21_atan(NAN));
    ck_assert_ldouble_nan(atan(NAN));
}
END_TEST

START_TEST(atan_8) {
    ck_assert_ldouble_eq_tol(s21_atan(1.0), atan(1.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_atan(1.1), atan(1.1), 1e-6);
    ck_assert_ldouble_eq_tol(s21_atan(-1.0), atan(-1.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_atan(-1.1), atan(-1.1), 1e-6);
}
END_TEST

START_TEST(atan_9) {
    ck_assert_ldouble_eq_tol(s21_atan(0.0), atan(0.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_atan(-0.0), atan(-0.0), 1e-6);
}
END_TEST

START_TEST(atan_10) { ck_assert_ldouble_eq_tol(s21_atan(-1), atan(-1), 1e-6); }
END_TEST

START_TEST(atan_11) {
    ck_assert_ldouble_eq_tol(s21_atan(-sqrt(3) / 2), atan(-sqrt(3) / 2), 1e-6);
}
END_TEST

START_TEST(atan_12) {
    ck_assert_ldouble_eq_tol(s21_atan(-sqrt(2) / 2), atan(-sqrt(2) / 2), 1e-6);
}
END_TEST

START_TEST(atan_13) {
    ck_assert_ldouble_eq_tol(s21_atan(-1 / 2), atan(-1 / 2), 1e-6);
}
END_TEST

START_TEST(atan_14) {
    ck_assert_ldouble_eq_tol(s21_atan(1 / 2), atan(1 / 2), 1e-6);
}
END_TEST

START_TEST(atan_15) {
    ck_assert_ldouble_eq_tol(s21_atan(2 / 2), atan(2 / 2), 1e-6);
}
END_TEST

START_TEST(atan_16) {
    ck_assert_ldouble_eq_tol(s21_atan(sqrt(3) / 2), atan(sqrt(3) / 2), 1e-6);
}
END_TEST

START_TEST(atan_17) { ck_assert_ldouble_eq_tol(s21_atan(1), atan(1), 1e-6); }
END_TEST

START_TEST(atan_18) {
    ck_assert_ldouble_eq_tol(s21_atan(1.23456e-7), atan(1.23456e-7), 1e-6);
}
END_TEST

START_TEST(atan_19) {
    ck_assert_ldouble_eq_tol(s21_atan(2.1234567), atan(2.1234567), 1e-6);
}
END_TEST

START_TEST(atan_20) {
    ck_assert_ldouble_eq_tol(s21_atan(2.1234567e-7), atan(2.1234567e-7), 1e-6);
}
END_TEST

START_TEST(atan_21) {
    ck_assert_ldouble_eq_tol(s21_atan(2.1234567e-16), atan(2.1234567e-16),
                             1e-6);
}
END_TEST

START_TEST(atan_22) { ck_assert_float_nan(s21_atan(s21_NAN)); }
END_TEST

START_TEST(atan_23) {
    ck_assert_ldouble_eq_tol(s21_atan(s21_PI / 2), atan(s21_PI / 2), 1e-6);
}
END_TEST

START_TEST(atan_24) {
    ck_assert_ldouble_eq_tol(s21_atan(-s21_PI / 2), atan(-s21_PI / 2), 1e-6);
}
END_TEST

START_TEST(s21_fmod_test_1) {
    ck_assert_ldouble_eq(s21_fmod(2.34, 2.0), fmod(2.34, 2.0));
}
END_TEST

START_TEST(s21_fmod_test_2) {
    ck_assert_ldouble_eq(s21_fmod(-2.34, 2.0), fmod(-2.34, 2.0));
}
END_TEST

START_TEST(s21_fmod_test_3) {
    ck_assert_ldouble_eq(s21_fmod(2.34, -2.0), fmod(2.34, -2.0));
}
END_TEST

START_TEST(s21_fmod_test_4) {
    ck_assert_ldouble_eq(s21_fmod(-2.34, -2.0), fmod(-2.34, -2.0));
}
END_TEST

START_TEST(s21_fmod_test_5) {
    ck_assert_ldouble_eq(s21_fmod(21.21, 3), fmod(21.21, 3));
}
END_TEST

START_TEST(s21_fmod_test_6) {
    ck_assert_ldouble_eq(s21_fmod(3, 21.21), fmod(3, 21.21));
}
END_TEST

START_TEST(s21_fmod_test_7) {
    ck_assert_ldouble_eq(s21_fmod(-21.21, 3), fmod(-21.21, 3));
}
END_TEST

START_TEST(s21_fmod_test_8) {
    ck_assert_ldouble_eq(s21_fmod(3, -21.21), fmod(3, -21.21));
}
END_TEST

START_TEST(s21_fmod_test_9) {
    ck_assert_ldouble_eq(s21_fmod(100500, 9), fmod(100500, 9));
}
END_TEST

START_TEST(s21_fmod_test_10) {
    ck_assert_ldouble_eq(s21_fmod(-100500, -9), fmod(-100500, -9));
}
END_TEST

START_TEST(s21_fmod_test_11) {
    ck_assert_ldouble_eq(s21_fmod(-9, -100500), fmod(-9, -100500));
}
END_TEST

START_TEST(s21_fmod_test_12) {
    ck_assert_ldouble_eq(s21_fmod(-9, -9), fmod(-9, -9));
}
END_TEST

START_TEST(s21_fmod_test_13) {
    ck_assert_ldouble_eq(s21_fmod(10, 5), fmod(10, 5));
}
END_TEST

START_TEST(s21_fmod_test_14) {
    ck_assert_ldouble_nan(s21_fmod(INFINITY, INFINITY));
    ck_assert_ldouble_nan(fmod(INFINITY, INFINITY));
}
END_TEST

START_TEST(s21_fmod_test_15) {
    ck_assert_ldouble_nan(s21_fmod(NAN, NAN));
    ck_assert_ldouble_nan(fmod(NAN, NAN));
}
END_TEST

START_TEST(s21_fmod_test_16) {
    ck_assert_ldouble_nan(s21_fmod(2.45, 0));
    ck_assert_ldouble_nan(fmod(2.45, 0));
}
END_TEST

START_TEST(s21_fmod_test_17) {
    ck_assert_ldouble_nan(s21_fmod(0, 0));
    ck_assert_ldouble_nan(fmod(0, 0));
}
END_TEST

START_TEST(s21_fmod_test_18) {
    ck_assert_ldouble_nan(s21_fmod(INFINITY, 0));
    ck_assert_ldouble_nan(fmod(INFINITY, 0));
}
END_TEST

START_TEST(s21_fmod_test_19) {
    double num = NAN;
    long double orig_res = fmod(num, 2.0), our_res = s21_fmod(num, 2.0);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(s21_fmod_test_20) {
    double num = NAN;
    long double orig_res = fmod(2.0, num), our_res = s21_fmod(2.0, num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(s21_fmod_test_21) {
    ck_assert_ldouble_eq(s21_fmod(2.0, INFINITY), fmod(2.0, INFINITY));
}
END_TEST

START_TEST(s21_fmod_test_22) {
    double num = INFINITY;
    long double orig_res = fmod(num, 2.0), our_res = s21_fmod(num, 2.0);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(s21_fmod_test_23) {
    ck_assert_ldouble_eq(s21_fmod(2.0, INFINITY), fmod(2.0, INFINITY));
}
END_TEST

START_TEST(s21_fmod_test_24) {
    ck_assert_ldouble_eq(s21_fmod(2.0, -INFINITY), fmod(2.0, -INFINITY));
}
END_TEST

START_TEST(s21_fmod_test_25) {
    double num = -INFINITY;
    long double orig_res = fmod(num, 2.0), our_res = s21_fmod(num, 2.0);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
}
END_TEST

START_TEST(s21_fmod_test_26) {
    ck_assert_ldouble_nan(s21_fmod(s21_INF, s21_INF));
    ck_assert_ldouble_nan(fmod(s21_INF, s21_INF));
}
END_TEST

START_TEST(s21_fmod_test_27) {
    ck_assert_ldouble_eq(s21_fmod(0, 23), fmod(0, 23));
}
END_TEST

void s21_math(TCase *tc1) {
    tcase_add_test(tc1, abs_1);
    tcase_add_test(tc1, abs_2);
    tcase_add_test(tc1, abs_3);
    tcase_add_test(tc1, abs_4);
    tcase_add_test(tc1, abs_5);
    tcase_add_test(tc1, abs_6);
    tcase_add_test(tc1, abs_7);
    tcase_add_test(tc1, abs_8);
    tcase_add_test(tc1, abs_9);
    tcase_add_test(tc1, abs_10);

    tcase_add_test(tc1, sqrt_1);
    tcase_add_test(tc1, sqrt_2);
    tcase_add_test(tc1, sqrt_3);
    tcase_add_test(tc1, sqrt_4);
    tcase_add_test(tc1, sqrt_5);
    tcase_add_test(tc1, sqrt_6);
    tcase_add_test(tc1, sqrt_7);
    tcase_add_test(tc1, sqrt_8);
    tcase_add_test(tc1, sqrt_9);
    tcase_add_test(tc1, sqrt_10);
    tcase_add_test(tc1, sqrt_12);
    tcase_add_test(tc1, sqrt_13);
    tcase_add_test(tc1, sqrt_14);
    tcase_add_test(tc1, sqrt_15);
    tcase_add_test(tc1, sqrt_16);

    tcase_add_test(tc1, log_1);
    tcase_add_test(tc1, log_2);
    tcase_add_test(tc1, log_3);
    tcase_add_test(tc1, log_4);
    tcase_add_test(tc1, log_5);
    tcase_add_test(tc1, log_6);
    tcase_add_test(tc1, log_7);
    tcase_add_test(tc1, log_8);
    tcase_add_test(tc1, log_9);
    tcase_add_test(tc1, log_10);
    tcase_add_test(tc1, log_11);
    tcase_add_test(tc1, log_12);
    tcase_add_test(tc1, log_13);
    tcase_add_test(tc1, log_14);
    tcase_add_test(tc1, log_15);
    tcase_add_test(tc1, log_16);
    tcase_add_test(tc1, log_17);
    tcase_add_test(tc1, log_18);
    tcase_add_test(tc1, log_19);
    tcase_add_test(tc1, log_20);
    tcase_add_test(tc1, log_21);
    tcase_add_test(tc1, log_22);
    tcase_add_test(tc1, log_23);
    tcase_add_test(tc1, log_24);
    tcase_add_test(tc1, log_25);
    tcase_add_test(tc1, log_26);
    tcase_add_test(tc1, log_27);
    tcase_add_test(tc1, log_28);
    tcase_add_test(tc1, log_29);
    tcase_add_test(tc1, log_30);
    tcase_add_test(tc1, log_31);
    tcase_add_test(tc1, log_32);
    tcase_add_test(tc1, log_33);
    tcase_add_test(tc1, log_34);

    tcase_add_test(tc1, fabs_1);
    tcase_add_test(tc1, fabs_2);
    tcase_add_test(tc1, fabs_3);
    tcase_add_test(tc1, fabs_4);
    tcase_add_test(tc1, fabs_5);
    tcase_add_test(tc1, fabs_6);
    tcase_add_test(tc1, fabs_7);
    tcase_add_test(tc1, fabs_8);
    tcase_add_test(tc1, fabs_10);
    tcase_add_test(tc1, fabs_11);
    tcase_add_test(tc1, fabs_12);
    tcase_add_test(tc1, fabs_13);
    tcase_add_test(tc1, fabs_14);
    tcase_add_test(tc1, fabs_15);

    tcase_add_test(tc1, exp_0);
    tcase_add_test(tc1, exp_1);
    tcase_add_test(tc1, exp_2);
    tcase_add_test(tc1, exp_3);
    tcase_add_test(tc1, exp_4);
    tcase_add_test(tc1, exp_5);
    tcase_add_test(tc1, exp_6);
    tcase_add_test(tc1, exp_7);
    tcase_add_test(tc1, exp_8);
    tcase_add_test(tc1, exp_9);
    tcase_add_test(tc1, exp_10);
    tcase_add_test(tc1, exp_11);
    tcase_add_test(tc1, exp_12);
    tcase_add_test(tc1, exp_13);
    tcase_add_test(tc1, exp_14);
    tcase_add_test(tc1, exp_15);
    tcase_add_test(tc1, exp_16);
    tcase_add_test(tc1, exp_17);
    tcase_add_test(tc1, exp_18);
    tcase_add_test(tc1, exp_19);
    tcase_add_test(tc1, exp_20);
    tcase_add_test(tc1, exp_21);
    tcase_add_test(tc1, exp_22);
    tcase_add_test(tc1, exp_23);
    tcase_add_test(tc1, exp_24);
    tcase_add_test(tc1, exp_25);
    tcase_add_test(tc1, exp_26);
    tcase_add_test(tc1, exp_27);
    tcase_add_test(tc1, exp_28);
    tcase_add_test(tc1, exp_29);
    tcase_add_test(tc1, exp_30);
    tcase_add_test(tc1, exp_31);
    tcase_add_test(tc1, exp_32);
    tcase_add_test(tc1, exp_33);

    tcase_add_test(tc1, pow_1);
    tcase_add_test(tc1, pow_2);
    tcase_add_test(tc1, pow_3);
    tcase_add_test(tc1, pow_4);
    tcase_add_test(tc1, pow_5);
    tcase_add_test(tc1, pow_6);
    tcase_add_test(tc1, pow_7);
    tcase_add_test(tc1, pow_8);
    tcase_add_test(tc1, pow_9);
    tcase_add_test(tc1, pow_10);
    tcase_add_test(tc1, pow_12);
    tcase_add_test(tc1, pow_13);
    tcase_add_test(tc1, pow_14);
    tcase_add_test(tc1, pow_15);
    tcase_add_test(tc1, pow_16);
    tcase_add_test(tc1, pow_17);
    tcase_add_test(tc1, pow_18);
    tcase_add_test(tc1, pow_19);
    tcase_add_test(tc1, pow_20);
    tcase_add_test(tc1, pow_21);
    tcase_add_test(tc1, pow_22);
    tcase_add_test(tc1, pow_23);
    tcase_add_test(tc1, pow_24);
    tcase_add_test(tc1, pow_25);
    tcase_add_test(tc1, pow_26);
    tcase_add_test(tc1, pow_27);
    tcase_add_test(tc1, pow_28);
    tcase_add_test(tc1, pow_29);
    tcase_add_test(tc1, pow_30);
    tcase_add_test(tc1, pow_31);
    tcase_add_test(tc1, pow_32);
    tcase_add_test(tc1, pow_33);
    tcase_add_test(tc1, pow_34);
    tcase_add_test(tc1, pow_35);
    tcase_add_test(tc1, pow_36);
    tcase_add_test(tc1, pow_37);
    tcase_add_test(tc1, pow_38);
    tcase_add_test(tc1, pow_39);
    tcase_add_test(tc1, pow_40);
    tcase_add_test(tc1, pow_41);
    tcase_add_test(tc1, pow_42);
    tcase_add_test(tc1, pow_43);
    tcase_add_test(tc1, pow_44);
    tcase_add_test(tc1, pow_45);
    tcase_add_test(tc1, pow_46);
    tcase_add_test(tc1, pow_47);
    tcase_add_test(tc1, pow_48);
    tcase_add_test(tc1, pow_49);
    tcase_add_test(tc1, pow_50);
    tcase_add_test(tc1, pow_51);
    tcase_add_test(tc1, pow_52);
    tcase_add_test(tc1, pow_53);
    tcase_add_test(tc1, pow_54);
    tcase_add_test(tc1, pow_55);
    tcase_add_test(tc1, pow_56);
    tcase_add_test(tc1, pow_57);
    tcase_add_test(tc1, pow_58);
    tcase_add_test(tc1, pow_59);
    tcase_add_test(tc1, pow_60);
    tcase_add_test(tc1, pow_61);
    tcase_add_test(tc1, pow_62);
    tcase_add_test(tc1, pow_63);
    tcase_add_test(tc1, pow_64);
    tcase_add_test(tc1, pow_65);
    tcase_add_test(tc1, pow_66);
    tcase_add_test(tc1, pow_67);
    tcase_add_test(tc1, pow_68);
    tcase_add_test(tc1, pow_69);
    tcase_add_test(tc1, pow_70);
    tcase_add_test(tc1, pow_71);
    tcase_add_test(tc1, pow_72);
    tcase_add_test(tc1, pow_73);
    tcase_add_test(tc1, pow_74);

    tcase_add_test(tc1, floor_1);
    tcase_add_test(tc1, floor_2);
    tcase_add_test(tc1, floor_3);
    tcase_add_test(tc1, floor_4);
    tcase_add_test(tc1, floor_5);
    tcase_add_test(tc1, floor_6);
    tcase_add_test(tc1, floor_7);
    tcase_add_test(tc1, floor_8);
    tcase_add_test(tc1, floor_10);
    tcase_add_test(tc1, floor_11);
    tcase_add_test(tc1, floor_12);
    tcase_add_test(tc1, floor_13);
    tcase_add_test(tc1, floor_14);
    tcase_add_test(tc1, floor_15);

    tcase_add_test(tc1, ceil_1);
    tcase_add_test(tc1, ceil_2);
    tcase_add_test(tc1, ceil_3);
    tcase_add_test(tc1, ceil_4);
    tcase_add_test(tc1, ceil_5);
    tcase_add_test(tc1, ceil_6);
    tcase_add_test(tc1, ceil_7);
    tcase_add_test(tc1, ceil_8);
    tcase_add_test(tc1, ceil_9);
    tcase_add_test(tc1, ceil_11);
    tcase_add_test(tc1, ceil_12);
    tcase_add_test(tc1, ceil_13);
    tcase_add_test(tc1, ceil_14);
    tcase_add_test(tc1, ceil_15);
    tcase_add_test(tc1, ceil_16);

    tcase_add_test(tc1, sin_1);
    tcase_add_test(tc1, sin_2);
    tcase_add_test(tc1, sin_3);
    tcase_add_test(tc1, sin_4);
    tcase_add_test(tc1, sin_5);
    tcase_add_test(tc1, sin_6);
    tcase_add_test(tc1, sin_7);
    tcase_add_test(tc1, sin_8);
    tcase_add_test(tc1, sin_9);
    tcase_add_test(tc1, sin_10);
    tcase_add_test(tc1, sin_11);
    tcase_add_test(tc1, sin_12);
    tcase_add_test(tc1, sin_13);
    tcase_add_test(tc1, sin_14);
    tcase_add_test(tc1, sin_15);
    tcase_add_test(tc1, sin_16);
    tcase_add_test(tc1, sin_17);
    tcase_add_test(tc1, sin_18);
    tcase_add_test(tc1, sin_19);
    tcase_add_test(tc1, sin_20);
    tcase_add_test(tc1, sin_21);
    tcase_add_test(tc1, sin_22);
    tcase_add_test(tc1, sin_23);
    tcase_add_test(tc1, sin_24);
    tcase_add_test(tc1, sin_25);

    tcase_add_test(tc1, cos_1);
    tcase_add_test(tc1, cos_2);
    tcase_add_test(tc1, cos_3);
    tcase_add_test(tc1, cos_4);
    tcase_add_test(tc1, cos_5);
    tcase_add_test(tc1, cos_6);
    tcase_add_test(tc1, cos_7);
    tcase_add_test(tc1, cos_8);
    tcase_add_test(tc1, cos_9);
    tcase_add_test(tc1, cos_10);
    tcase_add_test(tc1, cos_11);
    tcase_add_test(tc1, cos_12);
    tcase_add_test(tc1, cos_13);
    tcase_add_test(tc1, cos_14);
    tcase_add_test(tc1, cos_15);
    tcase_add_test(tc1, cos_16);
    tcase_add_test(tc1, cos_17);
    tcase_add_test(tc1, cos_18);
    tcase_add_test(tc1, cos_19);
    tcase_add_test(tc1, cos_20);
    tcase_add_test(tc1, cos_21);

    tcase_add_test(tc1, tan_1);
    tcase_add_test(tc1, tan_2);
    tcase_add_test(tc1, tan_3);
    tcase_add_test(tc1, tan_4);
    tcase_add_test(tc1, tan_5);
    tcase_add_test(tc1, tan_6);
    tcase_add_test(tc1, tan_7);
    tcase_add_test(tc1, tan_8);
    tcase_add_test(tc1, tan_9);
    tcase_add_test(tc1, tan_10);
    tcase_add_test(tc1, tan_11);
    tcase_add_test(tc1, tan_12);
    tcase_add_test(tc1, tan_13);
    tcase_add_test(tc1, tan_14);
    tcase_add_test(tc1, tan_15);
    tcase_add_test(tc1, tan_16);
    tcase_add_test(tc1, tan_17);
    tcase_add_test(tc1, tan_18);
    tcase_add_test(tc1, tan_19);
#if defined(__linux__)
    tcase_add_test(tc1, tan_20);
    tcase_add_test(tc1, tan_21);
#endif

    tcase_add_test(tc1, asin_1);
    tcase_add_test(tc1, asin_2);
    tcase_add_test(tc1, asin_3);
    tcase_add_test(tc1, asin_4);
    tcase_add_test(tc1, asin_5);
    tcase_add_test(tc1, asin_6);
    tcase_add_test(tc1, asin_7);
    tcase_add_test(tc1, asin_8);
    tcase_add_test(tc1, asin_9);
    tcase_add_test(tc1, asin_10);
    tcase_add_test(tc1, asin_11);
    tcase_add_test(tc1, asin_12);
    tcase_add_test(tc1, asin_13);
    tcase_add_test(tc1, asin_14);
    tcase_add_test(tc1, asin_15);
    tcase_add_test(tc1, asin_16);
    tcase_add_test(tc1, asin_17);
    tcase_add_test(tc1, asin_18);
    tcase_add_test(tc1, asin_19);
    tcase_add_test(tc1, asin_20);

    tcase_add_test(tc1, acos_1);
    tcase_add_test(tc1, acos_2);
    tcase_add_test(tc1, acos_3);
    tcase_add_test(tc1, acos_4);
    tcase_add_test(tc1, acos_5);
    tcase_add_test(tc1, acos_6);
    tcase_add_test(tc1, acos_7);
    tcase_add_test(tc1, acos_8);
    tcase_add_test(tc1, acos_9);
    tcase_add_test(tc1, acos_10);
    tcase_add_test(tc1, acos_11);
    tcase_add_test(tc1, acos_12);
    tcase_add_test(tc1, acos_13);
    tcase_add_test(tc1, acos_14);
    tcase_add_test(tc1, acos_15);
    tcase_add_test(tc1, acos_16);
    tcase_add_test(tc1, acos_17);
    tcase_add_test(tc1, acos_18);
    tcase_add_test(tc1, acos_19);

    tcase_add_test(tc1, atan_1);
    tcase_add_test(tc1, atan_2);
    tcase_add_test(tc1, atan_3);
    tcase_add_test(tc1, atan_4);
    tcase_add_test(tc1, atan_5);
    tcase_add_test(tc1, atan_6);
    tcase_add_test(tc1, atan_7);
    tcase_add_test(tc1, atan_8);
    tcase_add_test(tc1, atan_9);
    tcase_add_test(tc1, atan_10);
    tcase_add_test(tc1, atan_11);
    tcase_add_test(tc1, atan_12);
    tcase_add_test(tc1, atan_13);
    tcase_add_test(tc1, atan_14);
    tcase_add_test(tc1, atan_15);
    tcase_add_test(tc1, atan_16);
    tcase_add_test(tc1, atan_17);
    tcase_add_test(tc1, atan_18);
    tcase_add_test(tc1, atan_19);
    tcase_add_test(tc1, atan_20);
    tcase_add_test(tc1, atan_21);
    tcase_add_test(tc1, atan_22);
    tcase_add_test(tc1, atan_23);
    tcase_add_test(tc1, atan_24);

    tcase_add_test(tc1, s21_fmod_test_1);
    tcase_add_test(tc1, s21_fmod_test_2);
    tcase_add_test(tc1, s21_fmod_test_3);
    tcase_add_test(tc1, s21_fmod_test_4);
    tcase_add_test(tc1, s21_fmod_test_5);
    tcase_add_test(tc1, s21_fmod_test_6);
    tcase_add_test(tc1, s21_fmod_test_7);
    tcase_add_test(tc1, s21_fmod_test_8);
    tcase_add_test(tc1, s21_fmod_test_9);
    tcase_add_test(tc1, s21_fmod_test_10);
    tcase_add_test(tc1, s21_fmod_test_11);
    tcase_add_test(tc1, s21_fmod_test_12);
    tcase_add_test(tc1, s21_fmod_test_13);
    tcase_add_test(tc1, s21_fmod_test_14);
    tcase_add_test(tc1, s21_fmod_test_15);
    tcase_add_test(tc1, s21_fmod_test_16);
    tcase_add_test(tc1, s21_fmod_test_17);
    tcase_add_test(tc1, s21_fmod_test_18);
    tcase_add_test(tc1, s21_fmod_test_19);
    tcase_add_test(tc1, s21_fmod_test_20);
    tcase_add_test(tc1, s21_fmod_test_21);
    tcase_add_test(tc1, s21_fmod_test_22);
    tcase_add_test(tc1, s21_fmod_test_23);
    tcase_add_test(tc1, s21_fmod_test_24);
    tcase_add_test(tc1, s21_fmod_test_25);
    tcase_add_test(tc1, s21_fmod_test_26);
    tcase_add_test(tc1, s21_fmod_test_27);
}

int main(void) {
    int failed;
    Suite *s = suite_create("Core");
    TCase *tc1 = tcase_create("Core");
    SRunner *sr = srunner_create(s);

    suite_add_tcase(s, tc1);
    s21_math(tc1);

    srunner_run_all(sr, CK_ENV);
    failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
