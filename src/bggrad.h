void bggrad(uint8_t *img, int width, int height, uint8_t *out, uint8_t *bg, int step, int jump);
void bggrad_deshake(uint8_t *img, int width, int height, uint8_t *out, uint8_t *bg, int step, int jump, int deshake);
void bggrad_block_deshake(uint8_t *img, int width, int height, uint8_t *out, uint8_t *bg, int step, int jump, 
                          int deshake, int blocksize);
void bggrad_init(int _atan2bits, int _maxmag, double noise_variance);
void bggrad_noise(uint8_t *img, int width, int height, uint8_t *out, uint8_t *bg, int step, int jump,
                  uint8_t *prev, uint16_t *noise, int noise_scale, int noise_step);
void bggrad_noise_black(uint8_t *img, int width, int height, uint8_t *out, uint8_t *bg, int step, int jump,
                        uint8_t *prev, uint16_t *noise, int noise_scale, int noise_step, uint8_t *intensity);