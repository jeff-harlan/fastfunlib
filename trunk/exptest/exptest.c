/*
Test implementation of exponential function.

*/

#include <stdarg.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <sys/time.h>

// XXX -- really want a fixmul macro/function
#define SHIFT(x,prec) mpz_tdiv_q_2exp(x,x,prec)

#define MAX_SERIES_STEPS 8

mpz_t _exp_x;
mpz_t _exp_t;
mpz_t _exp_one;
mpz_t _exp_a;
mpz_t _exp_s0;
mpz_t _exp_s1;
mpz_t _exp_x2;

mpz_t _exp_pows[MAX_SERIES_STEPS];
mpz_t _exp_sums[MAX_SERIES_STEPS];

void fix_read_coeff_array(FILE *fp, mpz_t *coeffs, int *count);

double timing()
{
    double v;
    struct timeval t;
    gettimeofday(&t, NULL);
    v = (double) t.tv_usec;
    v = v + 1e6 * (double) t.tv_sec;
    return v;
}

void fix_init_data()
{
    int i;
    mpz_init(_exp_x);
    mpz_init(_exp_x2);
    mpz_init(_exp_one);
    mpz_init(_exp_t);
    mpz_init(_exp_a);
    mpz_init(_exp_s0);
    mpz_init(_exp_s1);
    for (i=0; i<MAX_SERIES_STEPS; i++)
    {
        mpz_init(_exp_pows[i]);
        mpz_init(_exp_sums[i]);
    }
}

void fix_clear_data()
{
    int i;
    mpz_clear(_exp_x);
    mpz_clear(_exp_x2);
    mpz_clear(_exp_one);
    mpz_clear(_exp_t);
    mpz_clear(_exp_a);
    mpz_clear(_exp_s0);
    mpz_clear(_exp_s1);
    for (i=0; i<MAX_SERIES_STEPS; i++)
    {
        mpz_clear(_exp_pows[i]);
        mpz_clear(_exp_sums[i]);
    }
}

void mpz_fixed_one(mpz_t x, int prec)
{
    mpz_set_ui(x, 1);
    mpz_mul_2exp(x, x, prec);
}

/*
First version -- uses Taylor series for exp(x) directly, with two steps.
*/
void fix_exp(mpz_t z, mpz_t x, int prec)
{
    int k;
    int r = 8;
    //prec += 20;
    //mpz_set(_exp_x, x);
    //mpz_tdiv_q_2exp(_exp_x, _exp_x, r);
    mpz_tdiv_q_2exp(_exp_x, x, r);
    mpz_set_ui(_exp_s0, 1);
    mpz_mul_2exp(_exp_s0, _exp_s0, prec);
    mpz_set(_exp_s1, _exp_s0);
    mpz_mul(_exp_x2, _exp_x, _exp_x);
    SHIFT(_exp_x2,prec);
    mpz_set(_exp_a, _exp_x2);
    k = 2;
    while(1)
    {
        mpz_tdiv_q_ui(_exp_a, _exp_a, k);
        if (mpz_sgn(_exp_a) == 0)
            break;
        mpz_add(_exp_s0, _exp_s0, _exp_a);
        k += 1;
        mpz_tdiv_q_ui(_exp_a, _exp_a, k);
        if (mpz_sgn(_exp_a) == 0)
            break;
        mpz_add(_exp_s1, _exp_s1, _exp_a);
        k += 1;
        mpz_mul(_exp_a, _exp_a, _exp_x2);
        SHIFT(_exp_a,prec);
        if (mpz_sgn(_exp_a) == 0)
            break;
    }
    mpz_mul(_exp_s1, _exp_s1, _exp_x);
    SHIFT(_exp_s1, prec);
    mpz_add(_exp_s0, _exp_s0, _exp_s1);
    for(k=0; k<r; k++)
    {
        mpz_mul(_exp_s0, _exp_s0, _exp_s0);
        SHIFT(_exp_s0, prec);
    }
    //mpz_tdiv_q_ui(z, _exp_s0, 20);
    mpz_set(z, _exp_s0);
}



/*
Computes the exponential / trigonometric series

  alt = 0  -- c = cosh(x), s = sinh(x)
  alt = 1  -- c = cos(x), s = sin(x)
  alt = 2  -- c = exp(x), s = n/a

using the cosh/sinh series. Parameters:

  prec -- 
  r    -- number of argument reductions
  J    -- number of partitions of the series

*/
void exp_series(mpz_t c, mpz_t s, mpz_t x, int prec, int r, int J, int alt)
{
    int i, k, wp;

    wp = prec + 2*r + 10;

    mpz_fixed_one(_exp_one, wp);

    /*   x / 2^r, adjusted to wp   */
    mpz_mul_2exp(_exp_x, x, wp-prec);
    mpz_tdiv_q_2exp(_exp_x, _exp_x, r);

    if (J < 1)
        J = 1;
    
    for (i=0; i<J; i++)
    {
        if (i == 0)
        {
            mpz_set(_exp_pows[i], _exp_one);
        }
        else if (i == 1)
        {
            mpz_mul(_exp_pows[i], _exp_x, _exp_x);
            mpz_tdiv_q_2exp(_exp_pows[i], _exp_pows[i], wp);
        }
        else
        {
            mpz_mul(_exp_pows[i], _exp_pows[i-1], _exp_pows[1]);
            mpz_tdiv_q_2exp(_exp_pows[i], _exp_pows[i], wp);
        }
        mpz_set_ui(_exp_sums[i], 0);
    }

    if (J == 1)
    {
        mpz_mul(_exp_x, _exp_x, _exp_x);
        mpz_tdiv_q_2exp(_exp_x, _exp_x, wp);
        mpz_set(_exp_a, _exp_x);
    }
    else
    {
        mpz_mul(_exp_x, _exp_pows[J-1], _exp_pows[1]);
        mpz_tdiv_q_2exp(_exp_x, _exp_x, wp);
        mpz_set(_exp_a, _exp_pows[1]);
    }

    k = 2;
    while (mpz_sgn(_exp_a) != 0)
    {
        for (i=0; i<J; i++)
        {
            mpz_tdiv_q_ui(_exp_a, _exp_a, (k-1)*k);
            if ((alt == 1) && (k & 2))
            {
                mpz_sub(_exp_sums[i], _exp_sums[i], _exp_a);
            }
            else
            {
                mpz_add(_exp_sums[i], _exp_sums[i], _exp_a);
            }
            k += 2;
        }
        mpz_mul(_exp_a, _exp_a, _exp_x);
        mpz_tdiv_q_2exp(_exp_a, _exp_a, wp);
    }

    for (i=1; i<J; i++)
    {
        mpz_mul(_exp_sums[i], _exp_sums[i], _exp_pows[i]);
        mpz_tdiv_q_2exp(_exp_sums[i], _exp_sums[i], wp);
    }

    mpz_set(c, _exp_one);
    for (i=0; i<J; i++)
    {
        mpz_add(c, c, _exp_sums[i]);
    }

    /*
    Repeatedly apply the duplication formula

      cosh(2*x) = 2*cosh(x)^2 - 1
      cos(2*x) = 2*cos(x)^2 - 1
      exp(2*x) = exp(x)^2
    */

    if (alt == 2)
    {
        /* s = sqrt(|1-c^2|) */
        mpz_mul_2exp(_exp_one, _exp_one, wp);
        mpz_mul(s, c, c);
        mpz_sub(s, _exp_one, s);
        mpz_abs(s, s);
        mpz_sqrt(s, s);
        mpz_add(c, c, s);
        for (i=0; i<r; i++)
        {
            mpz_mul(c, c, c);
            mpz_tdiv_q_2exp(c, c, wp);
        }
    }
    else
    {
        for (i=0; i<r; i++)
        {
            mpz_mul(c, c, c);
            mpz_tdiv_q_2exp(c, c, wp-1);
            mpz_sub(c, c, _exp_one);
        }
        /* s = sqrt(|1-c^2|) */
        mpz_mul_2exp(_exp_one, _exp_one, wp);
        mpz_mul(s, c, c);
        mpz_sub(s, _exp_one, s);
        mpz_abs(s, s);
        mpz_sqrt(s, s);
    }

    mpz_tdiv_q_2exp(c, c, wp-prec);
    mpz_tdiv_q_2exp(s, s, wp-prec);

}

void benchmark_optimize_exp()
{
    int REPS;
    int prec;
    int i, k, r, J, best_r, best_J;
    double best_time, best_here;
    double t1, t2, elapsed;
    double mpfr_time;
    int accuracy, min_accuracy;

    mpz_t x, y, dummy;
    mpfr_t mx, my;

    mpfr_init(mx);
    mpfr_init(my);

    mpz_init(x);
    mpz_init(y);
    mpz_init(dummy);

    printf(" prec   acc   J   r     mpfr     this   faster\n");

    for (prec=53; prec<4000; prec+=prec/4)
    {
        if (prec < 300)
            REPS = 50;
        else if (prec < 600)
            REPS = 20;
        else if (prec < 1200)
            REPS = 5;
        else
            REPS = 2;

        mpz_set_str(x, "3332663724254167", 10);
        mpz_mul_2exp(x, x, prec-53);
        //mpz_div_ui(x, x, 3);

        min_accuracy = prec;
        best_time = 1e100;
        best_r = 0;
        best_J = 0;

        mpfr_set_prec(mx, prec);
        mpfr_set_prec(my, prec);
        mpfr_set_d(mx, 0.37, GMP_RNDN);
        //mpfr_div_ui(mx, mx, 3, GMP_RNDN);

        mpfr_time = 100000.0;
        for (k=0; k<3; k++)
        {
            t1 = timing();
            for (k=0; k<REPS; k++)
            {
                mpfr_exp(my, mx, GMP_RNDN);
            }
            t2 = timing();
            elapsed = (t2-t1)/REPS;
            if (elapsed < mpfr_time)
                mpfr_time = elapsed;
        }

        for (J=0; J<MAX_SERIES_STEPS; J++)
        {
            for (r=0; r<prec/2; r++)
            {
                t1 = timing();
                for (k=0; k<REPS; k++)
                {
                    exp_series(y, dummy, x, prec, r, J, 2);
                }
                t2 = timing();
                elapsed = (t2-t1) / REPS;
                if (elapsed < best_time)
                {
                    best_time = elapsed;
                    best_r = r;
                    best_J = J;
                }

                mpfr_set_z(mx, y, GMP_RNDN);
                mpfr_div_2ui(mx, mx, prec, GMP_RNDN);
                //mpfr_printf("Value:  %Rf\n", mx);
                mpfr_sub(mx, mx, my, GMP_RNDN);
                mpfr_abs(mx, mx, GMP_RNDN);
                if (!mpfr_zero_p(mx))
                {
                    accuracy = -(int)mpfr_get_exp(mx)+1;
                    if (accuracy < min_accuracy)
                    {
                        min_accuracy = accuracy;
                    }
                }
            }
        }

        mpfr_time *= 1000;
        best_time *= 1000;

        printf("%5d %5d %3d %3d %8d %8d   %.3f\n", prec, min_accuracy, best_J, best_r,
            (int)mpfr_time, (int)best_time, mpfr_time/best_time);

    }

    mpfr_clear(mx);
    mpfr_clear(my);

    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(dummy);

}



int main()
{
    fix_init_data();

    benchmark_optimize_exp();

    fix_clear_data();
}

