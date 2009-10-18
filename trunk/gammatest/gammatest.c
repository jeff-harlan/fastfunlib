#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include <sys/time.h>
#include <math.h>

#define MAX_GAMMA_COEFF 3000
mpz_t g_rfac;
mpz_t g_one;
mpz_t gamma_coeff[MAX_GAMMA_COEFF];

mpz_t ta;
mpz_t tb;
mpz_t tc;
mpz_t td;

int gamma_coeff_prec = 0;
int gamma_max_coeff_index = 0;

void ffl_init()
{
    int k;
    mpz_init(g_rfac);
    mpz_init(g_one);
    mpz_init(ta);
    mpz_init(tb);
    mpz_init(tc);
    mpz_init(td);
    for (k=0; k<MAX_GAMMA_COEFF; k++)
    {
        mpz_init(gamma_coeff[k]);
    }
}

void ffl_clear()
{
    int k;
    mpz_clear(g_rfac);
    mpz_clear(g_one);
    mpz_clear(ta);
    mpz_clear(tb);
    mpz_clear(tc);
    mpz_clear(td);
    for (k=0; k<MAX_GAMMA_COEFF; k++)
    {
        mpz_clear(gamma_coeff[k]);
    }
}

void load_gamma_coefficients()
{
    FILE *fp;
    int k;
    int l;

    fp = fopen("gamma_data.txt", "rt");

    if ((fp == NULL) || (fscanf(fp, "%d", &gamma_coeff_prec) != 1))
    {
        printf("Could not open gamma_data.txt!\n");
        exit(1);
    }

    for (k=0; k<MAX_GAMMA_COEFF; k++)
    {
        if (gmp_fscanf(fp, "%Zx\n", gamma_coeff[k]) != 1)
            break;
    }

    gamma_max_coeff_index = k;

    fclose(fp);
}

void printx(char *s, mpz_t x, int prec)
{
    mpfr_t y;
    mpfr_init2(y, 53);
    mpfr_set_z(y, x, GMP_RNDN);
    mpfr_div_2ui(y, y, prec, GMP_RNDN);
    mpfr_printf("%s: %Rf\n", s, y);
    mpfr_clear(y);
}

/*
    gamma(x) = y * 2^n

    assumes x >= 0.5

    n is set to a nonzero value if x >> 10^0
*/
int gamma_taylor(mpz_t y, mpz_t x, int prec)
{
    int k, n, tmp, steps, terms, dprec;
    int wp;
    int expt = 0;
    wp = prec + 15;

    mpz_mul_2exp(ta, x, wp-prec);

    mpz_set_ui(g_one, 1);
    mpz_mul_2exp(g_one, g_one, wp);

    /* Reduce to [0.5,1.5) */
    mpz_tdiv_q_2exp(tb, ta, wp-1);
    n = mpz_get_si(tb);
    steps = (n-1)/2;
    if (steps)
    {
        /* TODO: using a polynomial of degree p,
           falling factorial can be evaluated using n/p+p
           full multiplications instead of n.
        */
        mpz_sub(ta, ta, g_one);
        mpz_set(g_rfac, ta);
        for (k=1; k<steps; k++)
        {
            mpz_sub(ta, ta, g_one);
            mpz_mul(g_rfac, g_rfac, ta);
            mpz_tdiv_q_2exp(g_rfac, g_rfac, wp);
            /* Don't grow too large */
            if (!(k % 4))
            {
                tmp = mpz_sizeinbase(g_rfac, 2) - wp;
                mpz_tdiv_q_2exp(g_rfac, g_rfac, tmp);
                expt += tmp;
            }
        }
    }
    else
    {
        mpz_set(g_rfac, g_one);
    }

    /* Polynomial is for G(1+x), so center on [-0.5,0.5) */
    mpz_sub(ta, ta, g_one);

    /* TODO: be both clever and correct here */
    if (wp < 1000)
    {
        terms = (int)(pow(wp, 0.76) + 2);
    }
    else
    {
        /* Valid up to at least 15000 bits */
        terms = (int)(pow(wp, 0.787) + 2);
    }

    dprec = gamma_coeff_prec-wp;
    mpz_tdiv_q_2exp(tb, gamma_coeff[terms], dprec);
    for (k=terms-1; k>=0; k--)
    {
        mpz_mul(tb, tb, ta);
        mpz_tdiv_q_2exp(tb, tb, wp);
        mpz_tdiv_q_2exp(tc, gamma_coeff[k], dprec);
        mpz_add(tb, tb, tc);
    }

    mpz_mul_2exp(g_rfac, g_rfac, wp - (wp-prec));
    mpz_div(y, g_rfac, tb);

    return expt;
}

double timing()
{
    double v;
    struct timeval t;
    gettimeofday(&t, NULL);
    v = (double) t.tv_usec;
    v = v + 1e6 * (double) t.tv_sec;
    return v;
}



void benchmark_gamma()
{
    int REPS;
    int prec;
    int i, k, r, J, best_r, best_J;
    double best_time, best_here;
    double t1, t2, elapsed;
    double mpfr_time;
    int accuracy, min_accuracy;
    int expt;

    mpz_t x, y, dummy;
    mpfr_t mx, my;

    mpfr_init(mx);
    mpfr_init(my);

    mpz_init(x);
    mpz_init(y);
    mpz_init(dummy);

    for (prec=53; prec<gamma_coeff_prec-100; prec+=prec/4)
    {
        if (prec < 300)
            REPS = 100;
        else if (prec < 600)
            REPS = 50;
        else if (prec < 1200)
            REPS = 10;
        else
            REPS = 2;

        mpz_set_ui(x, 57);
        mpz_mul_2exp(x, x, prec);
        mpz_div_ui(x, x, 10);

        min_accuracy = prec;
        best_time = 1e100;
        best_r = 0;
        best_J = 0;

        mpfr_set_prec(mx, prec+10);
        mpfr_set_prec(my, prec);
        mpfr_set_z(mx, x, GMP_RNDN);
        mpfr_div_2ui(mx, mx, prec, GMP_RNDN);

        mpfr_time = 1e100;
        for (i=0; i<10; i++)
        {
            t1 = timing();
            for (k=0; k<REPS; k++)
            {
                mpfr_gamma(my, mx, GMP_RNDN);
            }
            t2 = timing();
            elapsed = (t2-t1)/REPS;
            if (elapsed < mpfr_time)
                mpfr_time = elapsed;
        }

        best_time = 1e100;
        for (i=0; i<10; i++)
        {
            t1 = timing();
            for (k=0; k<REPS; k++)
            {
                expt = gamma_taylor(y, x, prec);
            }
            t2 = timing();
            elapsed = (t2-t1)/REPS;
            if (elapsed < best_time)
                best_time = elapsed;
        }

        mpfr_set_z(mx, y, GMP_RNDN);
        mpfr_mul_2ui(mx, mx, expt, GMP_RNDN);
        mpfr_div_2ui(mx, mx, prec, GMP_RNDN);

        mpfr_sub(mx, mx, my, GMP_RNDN);
        mpfr_div(mx, mx, my, GMP_RNDN);
        mpfr_abs(mx, mx, GMP_RNDN);

        if (!mpfr_zero_p(mx))
        {
            accuracy = -(int)mpfr_get_exp(mx)+1;
            if (accuracy < min_accuracy)
            {
                min_accuracy = accuracy;
            }
        }

        mpfr_time *= 1000;
        best_time *= 1000;
        printf("%5d %5d %10ld %10d   %.3f\n", prec, min_accuracy,
            (long)mpfr_time, (int)best_time, mpfr_time/best_time);

    }

    mpfr_clear(mx);
    mpfr_clear(my);

    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(dummy);

}



int main()
{
    ffl_init();
    load_gamma_coefficients();
    benchmark_gamma();
    ffl_clear();
}


