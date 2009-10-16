#include <stdarg.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <sys/time.h>

#define MAX_SERIES_STEPS 10

mpz_t _log_x;
mpz_t _log_t;
mpz_t _log_one;
mpz_t _log_a;
mpz_t _log_s0;
mpz_t _log_s1;
mpz_t _log_x2;

mpz_t _log_pows[MAX_SERIES_STEPS];
mpz_t _log_sums[MAX_SERIES_STEPS];

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
    mpz_init(_log_x);
    mpz_init(_log_x2);
    mpz_init(_log_one);
    mpz_init(_log_t);
    mpz_init(_log_a);
    mpz_init(_log_s0);
    mpz_init(_log_s1);
    for (i=0; i<MAX_SERIES_STEPS; i++)
    {
        mpz_init(_log_pows[i]);
        mpz_init(_log_sums[i]);
    }
}

void fix_clear_data()
{
    int i;
    mpz_clear(_log_x);
    mpz_clear(_log_x2);
    mpz_clear(_log_one);
    mpz_clear(_log_t);
    mpz_clear(_log_a);
    mpz_clear(_log_s0);
    mpz_clear(_log_s1);
    for (i=0; i<MAX_SERIES_STEPS; i++)
    {
        mpz_clear(_log_pows[i]);
        mpz_clear(_log_sums[i]);
    }
}

void mpz_fixed_one(mpz_t x, int prec)
{
    //mpz_set_ui(x, 1);
    //mpz_mul_2exp(x, x, prec);
    mpz_set_ui(x, 0);
    mpz_setbit(x, prec);
}

void log_series(mpz_t y, mpz_t x, int prec, int r, int J)
{
    int i, k, wp;

    wp = prec + r + 10;

    mpz_fixed_one(_log_one, wp);
    mpz_mul_2exp(_log_x, x, wp-prec);

    for (i=0; i<r; i++)
    {
        mpz_mul_2exp(_log_x, _log_x, wp);
        mpz_sqrt(_log_x, _log_x);
    }

    // x = (x-1)/(x+1)
    mpz_add(_log_t, _log_x, _log_one);
    mpz_sub(_log_x, _log_x, _log_one);
    mpz_mul_2exp(_log_x, _log_x, wp);
    mpz_tdiv_q(_log_x, _log_x, _log_t);

    if (J < 1)
        J = 1;

    for (i=0; i<J; i++)
    {
        if (i == 0)
        {
            mpz_set(_log_pows[i], _log_one);
        }
        else if (i == 1)
        {
            mpz_mul(_log_pows[i], _log_x, _log_x);
            mpz_tdiv_q_2exp(_log_pows[i], _log_pows[i], wp);
        }
        else
        {
            mpz_mul(_log_pows[i], _log_pows[i-1], _log_pows[1]);
            mpz_tdiv_q_2exp(_log_pows[i], _log_pows[i], wp);
        }
        mpz_set_ui(_log_sums[i], 0);
    }

    // Set a = x, _log_x = x^2
    if (J == 1)
    {
        mpz_set(_log_a, _log_x);
        mpz_mul(_log_x, _log_x, _log_x);
        mpz_tdiv_q_2exp(_log_x, _log_x, wp);
    }
    // Set a = x, _log_x = x^(2*J)
    else
    {
        mpz_set(_log_a, _log_x);
        mpz_mul(_log_x, _log_pows[J-1], _log_pows[1]);
        mpz_tdiv_q_2exp(_log_x, _log_x, wp);
    }

    // Main Taylor series loop
    k = 1;
    while (mpz_sgn(_log_a) != 0)
    {
        for (i=0; i<J; i++)
        {
            mpz_tdiv_q_ui(_log_t, _log_a, k);
            mpz_add(_log_sums[i], _log_sums[i], _log_t);
            k += 2;
        }
        mpz_mul(_log_a, _log_a, _log_x);
        mpz_tdiv_q_2exp(_log_a, _log_a, wp);
    }

    for (i=1; i<J; i++)
    {
        mpz_mul(_log_sums[i], _log_sums[i], _log_pows[i]);
        mpz_tdiv_q_2exp(_log_sums[i], _log_sums[i], wp);
    }

    mpz_set_ui(y, 0);
    for (i=0; i<J; i++)
    {
        mpz_add(y, y, _log_sums[i]);
    }

    mpz_tdiv_q_2exp(y, y, wp-prec-r-1);
}

void benchmark_optimize_log()
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

    for (prec=53; prec<6000; prec+=prec/4)
    {
        if (prec < 300)
            REPS = 100;
        else if (prec < 600)
            REPS = 50;
        else if (prec < 1200)
            REPS = 10;
        else
            REPS = 2;

        mpz_set_ui(x, 137);
        mpz_mul_2exp(x, x, prec);
        mpz_div_ui(x, x, 100);

        min_accuracy = prec;
        best_time = 1e100;
        best_r = 0;
        best_J = 0;

        mpfr_set_prec(mx, prec);
        mpfr_set_prec(my, prec);
        mpfr_set_str(mx, "1.37", 10, GMP_RNDN);

        mpfr_time = 1e100;
        for (i=0; i<10; i++)
        {
            t1 = timing();
            for (k=0; k<REPS; k++)
            {
                mpfr_log(my, mx, GMP_RNDN);
            }
            t2 = timing();
            elapsed = (t2-t1)/REPS;
            if (elapsed < mpfr_time)
                mpfr_time = elapsed;
        }

        for (J=1; J<MAX_SERIES_STEPS; J++)
        {
            for (r=3; r*r<prec+30; r++)
            {
                for (i=0; i<3; i++)
                {
                    t1 = timing();
                    for (k=0; k<REPS; k++)
                    {
                        log_series(y, x, prec, r, J);
                    }
                    t2 = timing();
                    elapsed = (t2-t1) / REPS;

                    if (elapsed < best_time)
                    {
                        best_time = elapsed;
                        best_r = r;
                        best_J = J;
                    }
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

    benchmark_optimize_log();

    fix_clear_data();
}

