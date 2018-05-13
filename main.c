#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "image.h"
#include <omp.h>

#define N 50000
#define RATIOX 1

int main(int argc, char *argv[]) {
    double start, finish;
    double duration;

    int test = 0;

    double complex E1[N] = {0}, E2[N] = {0}, E3[N] = {0},
            oE1[N] = {0}, oE2[N] = {0}, oE3[N] = {0},
            nE1[N] = {0}, nE2[N] = {0}, nE3[N] = {0};
    double complex C1, C2, C3;
    double k1, k2, k3, dk, dx = 1.0 / RATIOX, dr = 1, kk1 = 0.01;
    double f, x, r;
    double bright = 2, r1 = 5, r3 = 100;
    int nx, nr, cc, i;


//    omp_set_num_threads(atoi(argv[1]));

    // 设置线程数量
    int num_threads = 2;
    omp_set_num_threads(num_threads);


    k1 = 2;
    k2 = 2 * k1 + 0.005;
    k3 = 3 * k1 + 0.02;
    dk = k3 - k1 - k2;
    C1 = I * dx / (k1 * dr * dr);
    C2 = I * dx / (k2 * dr * dr);
    C3 = I * dx / (k3 * dr * dr);

    Image *image;
    image = image_new(1200, N);
    image_fill(image, 0);

    start = omp_get_wtime();

//启动并行域
#pragma omp parallel
    {

#pragma omp for private(r)
        for (nr = 1; nr <= N - 1; nr += 1) {
            r = nr * dr;
            oE1[nr] = exp(-pow((r - 150 * dr) / r1, 2));
            oE3[nr] = exp(-pow((r - 150 * dr) / r3, 2));
        }


#pragma omp for private(r)
        for (nr = 1; nr <= N - 2; nr += 1) {
            r = nr * dr;
            E1[nr] = oE1[nr] - C1 / 2 * (oE1[nr + 1] + oE1[nr - 1] - 2 * oE1[nr]);
            E3[nr] = oE3[nr] - C3 / 2 * (oE3[nr + 1] + oE3[nr - 1] - 2 * oE3[nr]);
        }



//所有线程达到同步
#pragma omp barrier

#pragma omp single
        nx = 1;

//        for (nx = 1; nx <= 2400 * RATIOX; nx++) {
//            for (nx = 1; nx <= 1000; nx++)
        while(nx <= 2400 * RATIOX) {

//单线程执行
#pragma omp single
            {
                test++;
                x = nx * dx;
            }

//并行执行for循环，有隐性的barrier
#pragma omp for private(r, f)
            for (nr = 1; nr <= N - 2; nr++) {
                r = nr * dr;
                f = 0;
                if (nx >= 600 * RATIOX && nx <= 610 * RATIOX) f = 30;
                nE1[nr] = oE1[nr] - C1 * (E1[nr + 1] + E1[nr - 1] - 2 * E1[nr]);
                nE2[nr] = oE2[nr] - C2 * (E2[nr + 1] + E2[nr - 1] - 2 * E2[nr]) -
                          0.5 * I * dx * kk1 * f * E3[nr] * conj(E1[nr]) * cexp(-I * dk * x);
                nE3[nr] = oE3[nr] - C3 * (E3[nr + 1] + E3[nr - 1] - 2 * E3[nr]);

            cc = cabs(nE1[nr]) * 200;
            if (cc > 255) cc = 255;
            if (cc < 0) cc = 0;
            if (nr % 2 == 0 && (nx / 2) % RATIOX == 0) image_set_pixel(image, nx / 2.0 / RATIOX, nr / 3.0, cc);
            cc = cabs(nE2[nr]) * 100 * bright;
            if (cc > 255) cc = 255;
            if (cc < 0) cc = 0;
            if (nr % 2 == 0 && (nx / 2) % RATIOX == 0)
                image_set_pixel(image, nx / 2.0 / RATIOX, nr / 3.0 + N / 3, cc);
            cc = cabs(nE3[nr]) * 100 * bright;
            if (cc > 255) cc = 255;
            if (cc < 0) cc = 0;
            if (nr % 2 == 0 && (nx / 2) % RATIOX == 0)
                image_set_pixel(image, nx / 2.0 / RATIOX, nr / 3.0 + N / 3 * 2, cc * (256 + 1));
            }


//单线程执行
#pragma omp single
            {
                nE1[0] = nE1[1] * 0;
                nE1[N - 1] = 0;
                nE2[0] = 0;
                nE2[N - 1] = 0;
                nE3[0] = 0;
                nE3[N - 1] = 0;
            }

//并行执行for循环，有隐性的barrier
//#pragma omp single
//            {
#pragma omp for
            for (i = 0; i < N; i++) {
                    oE1[i] = E1[i];
                    oE2[i] = E2[i];
                    oE3[i] = E3[i];
                    E1[i] = nE1[i];
                    E2[i] = nE2[i];
                    E3[i] = nE3[i];
                }
#pragma omp single
            nx++;
        }

//并行域结束
    }


    //结束时间
    finish = omp_get_wtime();

    //实际运行时间
    duration = finish - start;

    image_save(image, "pic.pgm");

    image_free(image);

    printf("Number of threads in environment: %d\n", omp_get_max_threads());
    printf("Duration: %f seconds\n", duration);
    printf("Speed: %f GDots/s\n", 2400 * RATIOX * 3 * N * 1.0e-009/ duration);



    return 0;

}
