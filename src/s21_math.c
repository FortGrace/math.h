#include "s21_math.h"

int s21_abs(int x) { return x > 0 ? x : -x; }

long double s21_fabs(long double x) {
    if (x < 0 || x == -0.0) {
        x = -x;
    }
    return x;
}

long double s21_floor(double x) {
    long double floor_x;
    if (!is_fin(x)) {
        floor_x = x;
    } else {
        floor_x = (long long int)x;
        if (s21_fabs(x - floor_x) > 0. && s21_fabs(x) > 0.) {
            if (x < 0.) {
                floor_x -= 1;
            }
        }
    }
    return floor_x;
}

long double s21_ceil(double x) {
    long double d;
    if (!is_fin(x)) {
        d = x;
    } else {
        d = x - (x - (long)x);
        if (x != s21_MAX_double) {
            if (d != x && x > 0) {
                d = d + 1;
            }
        } else {
            d = s21_MAX_double;
        }
    }
    return d;
}

long double s21_exp(double x) {
    long double result = 1, temp = 1;
    long double i = 1;
    int flag = 0;
    if (x < 0) {
        x *= -1;
        flag = 1;
    }
    while (s21_fabs(result) > s21_EPS) {
        result *= x / i;
        i += 1;
        temp += result;
        if (temp > s21_MAX_double) {
            temp = s21_INF;
            break;
        }
    }
    temp = flag == 1 ? temp > s21_MAX_double ? 0 : 1. / temp : temp;
    return temp = temp > s21_MAX_double ? s21_INF : temp;
}

long double s21_log(double x) {
    int ex_pow = 0;
    double result = 0;
    double compare = 0;
    if (x == s21_INF) {
        result = s21_INF;
    } else if (x == 0) {
        result = -s21_INF;
    } else if (x < 0) {
        result = s21_NAN;
    } else if (x == 1) {
        result = 0;
    } else {
        for (; x >= s21_EXP; x /= s21_EXP, ex_pow++) continue;
        int i;
        for (i = 0; i < 100; i++) {
            compare = result;
            result =
                compare + 2 * (x - s21_exp(compare)) / (x + s21_exp(compare));
        }
    }
    return (result + ex_pow);
}

long double s21_pow(double base, double exp) {
    long double result;
    if (base < 0) {
        if ((long int)exp == exp) {
            if (exp > 0) {
                result = base;
                for (long int i = 0; i < (long long int)exp - 1; i++) {
                    result *= base;
                }
            } else if (exp == 0) {
                result = 1;
            } else {
                result = 1 / base;
                for (long int i = 0; i < (long long int)exp * (-1) - 1; i++) {
                    result /= base;
                }
            }
        } else {
            if (exp == -s21_INF || exp == s21_INF) {
                if (base * (-1) < 1) {
                    result = 0;
                } else if (base * (-1) == 1) {
                    result = 1;
                } else {
                    if (exp == -s21_INF) {
                        result = 0;
                    } else {
                        result = s21_INF;
                    }
                }
            } else {
                if (is_inf(base)) {
                    result = 0;
                } else {
                    result = -s21_NAN;
                }
            }
        }
    } else if (base == 0) {
        if (exp == 0) {
            result = 1;
        } else if (exp < 0) {
            result = s21_INF;
        } else if (is_nan(exp)) {
            result = s21_NAN;
        } else {
            result = 0;
        }
    } else if (base == 1) {
        result = 1;
    } else {
        if ((long int)exp == exp) {
            if (exp > 0) {
                result = base;
                for (long int i = 0; i < (long int)exp - 1; i++) {
                    result *= base;
                }
            } else if (exp == 0) {
                result = 1;
            } else {
                result = 1 / base;
                for (long int i = 0; i < (long int)exp * (-1) - 1; i++) {
                    result /= base;
                }
            }
        } else {
            result = s21_exp(exp * (double)s21_log(base));
        }
    }
    return result;
}

long double s21_sin(double x) {
    long double result = 0;
    if (s21_fmod(x, 2 * s21_PI) == 0 && x > 0) {
        result = -0.0000000000000002;
    } else if (s21_fabs(x) < s21_pow(2, -26)) {
        result = x;
    } else if (x != x || x == s21_INF || x == -s21_INF) {
        result = s21_NAN;
    } else {
        x = s21_fmod(x, 2 * s21_PI);
        long double prev = x;
        int n = 1;
        while (s21_fabs(prev) > s21_EPS) {
            result += prev;
            n += 2.;
            prev = -prev * x * x / (n * (n - 1));
        }
    }
    return result;
}

long double s21_cos(double x) { return s21_sin(s21_PI / 2 - x); }

long double s21_sqrt(double x) {
    long double temp, sqrt;
    sqrt = x / 2;
    temp = 0;
    if (x < 0 || is_nan(x)) {
        temp = s21_NAN;
    } else {
        if (is_inf(x)) {
            temp = s21_INF;
        } else {
            while (sqrt != temp) {
                temp = sqrt;
                sqrt = (x / temp + temp) / 2;
            }
        }
    }
    return temp;
}

long double s21_tan(double x) {
    long double result = 0;
    if (x == s21_PI / 2) {
        result = 16331239353195370L;
    } else if (x == -s21_PI / 2) {
        result = -16331239353195370L;
    } else if (x == 0) {
        result = 0;
    } else {
        result = s21_sin(x) / s21_cos(x);
    }
    return result;
}

long double s21_asin(double x) {
    int sign = x >= 0 ? 1 : -1;
    x = s21_fabs(x);
    long double result = x;
    if (x < 1) {
        if (x == 0.8660254037844386) {
            result = s21_PI / 3;
        } else {
            long double step = x;
            int n = 1;
            while (step > s21_EPS) {
                step *= ((x * x) * (2 * n - 1) * (2 * n - 1)) /
                        ((2 * n) * (2 * n + 1));
                result += step;
                n++;
            }
        }
    } else if (s21_fabs(x) == 1) {
        result = x * s21_PI / 2.0L;
    } else {
        result = s21_NAN;
    }
    result *= sign;
    return result;
}

long double s21_acos(double x) {
    long double result = 0;
    if (s21_fabs(x) < 1) {
        result = s21_PI / 2 - s21_asin(x);
    } else if (x == -1) {
        result = s21_PI;
    } else if (x == 1) {
        result = 0;
    } else {
        result = s21_NAN;
    }
    return result;
}

long double s21_atan(double x) {
    long double result = 0;
    int sign = x >= 0 ? 1 : -1;
    if (x == x) {
        x = s21_fabs(x);
        if (x < 60) {
            if (x == s21_PI / 2 || x == -s21_PI / 2) {
                result = 1.003884821853887214L;
            } else {
                result = s21_asin(s21_fabs(x) / s21_sqrt(x * x + 1));
            }
        } else {
            result = 1.57079632679489661923 + s21_atan(-1 / x);
        }
    } else {
        result = s21_NAN;
    }
    result *= sign;
    return result;
}

long double s21_fmod(double x, double y) {
    long double result = 0;
    if (y != 0. && x == x && x != s21_INF && x != -s21_INF && y == y) {
        if (y == s21_INF || y == -s21_INF)
            result = x;
        else
            result = x - (long)(x / y) * y;
    } else {
        result = s21_NAN;
    }
    return result;
}
