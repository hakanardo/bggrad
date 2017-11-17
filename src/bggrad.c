#include <stdint.h>
#include "bggrad.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <alloca.h>

#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))

static uint8_t *atan2tab;
static int atan2bits;
static uint8_t *pdftab;
static int maxmag;

double phi(double x) {
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

double apdf(double a, double SNR2) {
    return exp(-(1.0/2.0)*SNR2)*(sqrt(M_PI) + sqrt(2.0*SNR2)*cos(a)*phi(sqrt(SNR2)*cos(a))*
                                              exp((1.0/2.0)*SNR2*pow(cos(a),2.0))*M_PI)/pow(2.0*M_PI,3.0/2.0);
}

void bggrad_init(int _atan2bits, int _maxmag, double noise_variance) {
    atan2bits = _atan2bits;
    maxmag = _maxmag;

    int atan2size = 1 << atan2bits;
    int offset = atan2size >> 1;
    int x, y;
    atan2tab = malloc(atan2size * atan2size);
    for (y=0; y<atan2size; y++) {
        for (x=0; x<atan2size; x++) {
            atan2tab[ (y << atan2bits) + x ] = atan2(y - offset, x - offset) * 256.0 / (2.0 * M_PI);
        }
    }

    int mag, grad;
    double pbg, pfg;
    pdftab = malloc((maxmag+1) * 256);
    pfg = 1.0 / (2.0 * M_PI);
    for (mag=0; mag<=maxmag; mag++) {
        for (grad=0; grad<256; grad++) {
            pbg = apdf( ((double) grad) * 2*M_PI / 256.0, ((double) mag) / noise_variance );
            pdftab[mag * 256 + grad] = pfg / (pfg + pbg) * 255.0;
        }
    }

}

void bggrad(uint8_t *img, int width, int height, uint8_t *out, uint8_t *bg, int step, int jump) {
    int x, y;
    int upcnt = rand() % jump;

    for (y=1; y<height-1; y++) {
        for (x=1; x<width-1; x++) {
            int dx = img[y*width + x - 1] - img[y*width + x + 1];
            int dy = img[(y-1)*width + x] - img[(y+1)*width + x];
            int mag = dx*dx + dy*dy;

            int xi = (dx + 0x100) >> (9 - atan2bits);
            int yi = (dy + 0x100) >> (9 - atan2bits);
            uint8_t grad = atan2tab[ (yi << atan2bits) + xi ];

            uint8_t diff = (grad - bg[y*width + x]) & 0xFF; 
            out[y*width + x] = pdftab[min(mag, maxmag) * 256 + diff];

            upcnt++;
            if (upcnt == jump) {
                upcnt = 0;
                if (diff >= 0x80) {
                    bg[y*width + x] -= step;
                } else if (diff > 0) {
                    bg[y*width + x] += step;
                }
            }
        }
    }
}

void bggrad_deshake(uint8_t *img, int width, int height, uint8_t *out, uint8_t *bg, int step, int jump, 
                    int deshake) {
    int x, y, dx, dy;
    int upcnt = rand() % jump;

    assert(deshake >= 1);

    for (y=deshake; y<height-deshake; y++) {
        for (x=deshake; x<width-deshake; x++) {
            int dx = img[y*width + x - 1] - img[y*width + x + 1];
            int dy = img[(y-1)*width + x] - img[(y+1)*width + x];
            int mag = dx*dx + dy*dy;

            int xi = (dx + 0x100) >> (9 - atan2bits);
            int yi = (dy + 0x100) >> (9 - atan2bits);
            uint8_t grad = atan2tab[ (yi << atan2bits) + xi ];

            uint8_t diff, mindiff=255;
            for (dy=-deshake; dy<=deshake; dy+=1) {
                for (dx=-deshake; dx<=deshake; dx+=1) {
                    diff = abs(grad - bg[(y+dy)*width + (x+dx)]) & 0xFF;
                    mindiff = min(diff, mindiff);
                }
            }
            out[y*width + x] = pdftab[min(mag, maxmag) * 256 + mindiff];

            upcnt++;
            if (upcnt == jump) {
                upcnt = 0;
                diff = (grad - bg[(y)*width + (x)]) & 0xFF; 
                if (diff >= 0x80) {
                    bg[y*width + x] -= step;
                } else if (diff > 0) {
                    bg[y*width + x] += step;
                }
            }
        }
    }
}

void bggrad_block_deshake(uint8_t *img, int width, int height, uint8_t *out, uint8_t *bg, int step, int jump, 
                          int deshake, int blocksize) {
    int x0, y0, x, y, dx, dy;
    int upcnt = rand() % jump;
    int *mag = alloca(width*height*sizeof(int));
    uint8_t *grad = alloca(width*height);

    for (y=1; y<height-1; y++) {
        for (x=1; x<width-1; x++) {
            int dx = img[y*width + x - 1] - img[y*width + x + 1];
            int dy = img[(y-1)*width + x] - img[(y+1)*width + x];
            mag[y*width + x] = dx*dx + dy*dy;

            int xi = (dx + 0x100) >> (9 - atan2bits);
            int yi = (dy + 0x100) >> (9 - atan2bits);
            uint8_t xygrad = grad[y*width + x] = atan2tab[ (yi << atan2bits) + xi ];

            
            if (step > 1) {
                upcnt++;
                if (upcnt == jump) {
                    upcnt = 0;
                    uint8_t diff = (xygrad - bg[(y)*width + (x)]) & 0xFF; 
                    if (diff >= 0x80) {
                        bg[y*width + x] -= step;
                    } else if (diff > 0) {
                        bg[y*width + x] += step;
                    }
                }
            }
        }
    }

    int dxsum=0, dysum=0, nsum=0;
    for (y0=deshake; y0<height-deshake-blocksize; y0+=blocksize) {
        for (x0=deshake; x0<width-deshake-blocksize; x0+=blocksize) {
            int bestdiff = 1<<30, bestdx, bestdy;
            for (dy=-deshake; dy<=deshake; dy+=1) {
                for (dx=-deshake; dx<=deshake; dx+=1) {
                    int diff = 0;
                    for (y=y0; y<y0+blocksize; y++) {
                        for (x=x0; x<x0+blocksize; x++) {
                            diff += abs(grad[y*width + x] - bg[(y+dy)*width + (x+dx)]);
                        }
                    }
                    if (diff < bestdiff) {
                        bestdiff = diff;
                        bestdx = dx;
                        bestdy = dy;
                    }
                }
            }
            dxsum += bestdx;
            dysum += bestdy;
            nsum += 1;

            for (y=y0; y<y0+blocksize; y++) {
                for (x=x0; x<x0+blocksize; x++) {
                    uint8_t diff = (grad[y*width + x] - bg[(y+bestdy)*width + (x+bestdx)]) & 0xFF; 
                    out[y*width + x] = pdftab[min(mag[y*width + x], maxmag) * 256 + diff];
                    if (step <= 1) {
                        upcnt++;
                        if (upcnt == jump) {
                            upcnt = 0;
                            if (diff >= 0x80) {
                                bg[(y+bestdy)*width + (x+bestdx)] -= step;
                            } else if (diff > 0) {
                                bg[(y+bestdy)*width + (x+bestdx)] += step;
                            }
                        }
                    }
                }
            }
        }
    }
    //printf("%d, %d\n", dxsum/nsum, dysum/nsum);
}

void bggrad_noise(uint8_t *img, int width, int height, uint8_t *out, uint8_t *bg, int step, int jump,
                  uint8_t *prev, uint16_t *noise, int noise_scale, int noise_step) {
    int i, x, y;
    int upcnt = rand() % jump;
    unsigned long scale = noise_scale * ((int) (256.0 * 256.0 / 1.48260221851 / 1.48260221851 * 4));

    for (y=1; y<height-1; y++) {
        for (x=1; x<width-1; x++) {
            int dx = img[y*width + x - 1] - img[y*width + x + 1];
            int dy = img[(y-1)*width + x] - img[(y+1)*width + x];
            unsigned long mag = dx*dx + dy*dy;

            int xi = (dx + 0x100) >> (9 - atan2bits);
            int yi = (dy + 0x100) >> (9 - atan2bits);
            uint8_t grad = atan2tab[ (yi << atan2bits) + xi ];

            uint8_t diff = (grad - bg[y*width + x]) & 0xFF; 
            
            int q1 = noise[y*width + x - 1] + 1*256;
            int q2 = noise[y*width + x + 1] + 1*256;
            int q3 = noise[(y-1)*width + x] + 1*256;
            int q4 = noise[(y+1)*width + x] + 1*256;
            unsigned long s2 = q1*q1 + q2*q2 + q3*q3 + q4*q4;
            unsigned long snr = (scale) * mag / (s2);
            //assert(fabs(snr - ((double) scale) * ((double) mag) / ((double) s2)) < 2);
            assert(snr >= 0);
            out[y*width + x] = pdftab[min(snr , maxmag) * 256 + diff];

            /*
            if (x==50 && y==25) {
                img[y*width + x] = 255;
                printf("%d\n", grad);
            }
            */

            upcnt++;
            if (upcnt == jump) {
                upcnt = 0;
                if (diff >= 0x80) {
                    bg[y*width + x] -= step;
                } else if (diff > 0) {
                    bg[y*width + x] += step;
                }
            }

        }
    }
    if (prev) {
        for (i=0; i<width*height; i++) {            
            int d = abs(img[i] - prev[i]) * 256;
            if (d > noise[i]) {
                noise[i] += noise_step;
            } else if (d < noise[i]) {
                noise[i] = max(noise[i] - noise_step, 2);
            }
        }
    }
}

void bggrad_noise_black(uint8_t *img, int width, int height, uint8_t *out, uint8_t *bg, int step, int jump,
                        uint8_t *prev, uint16_t *noise, int noise_scale, int noise_step, uint8_t *intensity) {
    int i, x, y;
    int upcnt = rand() % jump;
    unsigned long scale = noise_scale * ((int) (256.0 * 256.0 / 1.48260221851 / 1.48260221851 * 4));

    for (y=1; y<height-1; y++) {
        for (x=1; x<width-1; x++) {
            int dx = img[y*width + x - 1] - img[y*width + x + 1];
            int dy = img[(y-1)*width + x] - img[(y+1)*width + x];
            unsigned long mag = dx*dx + dy*dy;

            int xi = (dx + 0x100) >> (9 - atan2bits);
            int yi = (dy + 0x100) >> (9 - atan2bits);
            uint8_t grad = atan2tab[ (yi << atan2bits) + xi ];
            //if (x==54 && y==31) printf("%d\n", grad);

            uint8_t diff = (grad - bg[y*width + x]) & 0xFF; 
            
            int q1 = noise[y*width + x - 1] + 1*256;
            int q2 = noise[y*width + x + 1] + 1*256;
            int q3 = noise[(y-1)*width + x] + 1*256;
            int q4 = noise[(y+1)*width + x] + 1*256;
            unsigned long s2 = q1*q1 + q2*q2 + q3*q3 + q4*q4;
            unsigned long snr = (scale) * mag / (s2);
            //assert(fabs(snr - ((double) scale) * ((double) mag) / ((double) s2)) < 2);
            assert(snr >= 0);
            out[y*width + x] = pdftab[min(snr , maxmag) * 256 + diff];

            /*
            int minint = min(min(intensity[y*width + x - 1], intensity[y*width + x + 1]),
                             min(intensity[(y-1)*width + x], intensity[(y+1)*width + x]));
            int th = (q1 + q2 + q3 + q4) * 4 / 4 / 256;
            if (minint < th) out[y*width + x] = 127;
            */
            int x1 = min(img[y*width + x - 1], 255-img[y*width + x - 1]);
            int x2 = min(img[y*width + x + 1], 255-img[y*width + x + 1]);
            int y1 = min(img[(y-1)*width + x], 255-img[(y-1)*width + x]);
            int y2 = min(img[(y+1)*width + x], 255-img[(y+1)*width + x]);
            int mmint = min(max(x1, x2), max(y1, y2));
            int th = (q1 + q2 + q3 + q4) * 2 * 1.48260221851 / 4 / 256;
            if (mmint > th) {
                upcnt++;
                if (upcnt == jump) {
                    upcnt = 0;
                    if (diff >= 0x80) {
                        bg[y*width + x] -= step;
                    } else if (diff > 0) {
                        bg[y*width + x] += step;
                    }
                }
            }

        }
    }
    if (prev) {
        jump *= 4;
        upcnt = rand() % jump;
        for (i=0; i<width*height; i++) {            
            int d = abs(img[i] - prev[i]) * 256;
            if (d > noise[i]) {
                noise[i] += noise_step;
            } else if (d < noise[i]) {
                noise[i] = max(noise[i] - noise_step, 2);
            }
            /*
            upcnt++;
            if (upcnt == jump) {
                upcnt = 0;
                if (img[i] > intensity[i]) {
                    intensity[i] += step;
                } else if (img[i] < intensity[i]) {
                    intensity[i] -= step;
                }
            }
            */
        }
    }
    //i = 31*width + 136-1; printf("%d +- %f\n", intensity[i], noise[i]/256.0);
}