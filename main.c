#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "image.h"

const static int N = 768;
const static int RATIOX = 4;

int main() {
    double complex E1[N] = {0}, E2[N] = {0}, E3[N] = {0},
            oE1[N] = {0}, oE2[N] = {0}, oE3[N] = {0},
            nE1[N] = {0}, nE2[N] = {0}, nE3[N] = {0};
    double complex C1, C2, C3;
    double complex t;
    double k1, k2, k3, dk, dx = 1.0 / RATIOX, dr = 1, kk1 = 0.01;
    double f, x, r;
    double bright = 2, r1 = 5, r3 = 100;
    int nx, nr, cc;
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

    for (nr = 0; nr <= N - 1; nr++) {
        r = nr * dr;
        oE1[nr] = exp(-pow((r - 150 * dr) / r1, 2));
        oE3[nr] = exp(-pow((r - 150 * dr) / r3, 2));
    }
    for (nr = 1; nr <= N - 2; nr++) {
        r = nr * dr;
        E1[nr] = oE1[nr] - C1 / 2 * (oE1[nr + 1] + oE1[nr - 1] - 2 * oE1[nr]);
        E3[nr] = oE3[nr] - C3 / 2 * (oE3[nr + 1] + oE3[nr - 1] - 2 * oE3[nr]);
    }
    for (nx = 1; nx <= 2400 * RATIOX; nx++) {
        x = nx * dx;
        for (nr = 1; nr <= N - 2; nr++) {
            r = nr * dr;
            f = 0;
            if (nx >= 600 * RATIOX && nx <= 610 * RATIOX) f = 30;
            nE1[nr] = oE1[nr] - C1 * (E1[nr + 1] + E1[nr - 1] - 2 * E1[nr]);
            //nE2[nr] has a problem
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
        nE1[0] = nE1[1] * 0;
        nE1[N - 1] = 0;
        nE2[0] = 0;
        nE2[N - 1] = 0;
        nE3[0] = 0;
        nE3[N - 1] = 0;

        memcpy(oE1, E1, sizeof(
                double complex)*N);
        memcpy(oE2, E2, sizeof(
                double complex)*N);
        memcpy(oE3, E3, sizeof(
                double complex)*N);
        memcpy(E1, nE1, sizeof(
                double complex)*N);
        memcpy(E2, nE2, sizeof(
                double complex)*N);
        memcpy(E3, nE3, sizeof(
                double complex)*N);

    }

    image_save(image, "pic.pgm");

    image_free(image);


//    Image *image_fun;
//
//    image_fun = image_new (800, 800);
//
//    image_fill (image_fun, 0xaa);
//    draw_Taijitu (image_fun, 300, 0);
//    image_save (image_fun, "taiji_6.pgm");
//
//    image_free (image_fun);

    return 0;


}


