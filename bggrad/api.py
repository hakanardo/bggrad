from __future__ import division
from bggrad._bggrad import lib, ffi
import numpy as np

class ESTIMATE:
    pass

class BgGrad(object):
    noise_variance = None
    noise_scale = 10

    def __init__(self, maxfcnt=10000, deshake=0, blocksize=1, noise_variance=ESTIMATE, exposure_check=False):
        self.maxfcnt = maxfcnt
        self.bg = None
        self.fcnt = 0
        assert deshake >= 0
        self.deshake = deshake
        assert blocksize >= 1
        self.blocksize = blocksize
        self.exposure_check = exposure_check
        if self.noise_variance is None:
            if noise_variance is ESTIMATE:
                lib.bggrad_init(8, 10000, self.noise_scale)
            else:
                lib.bggrad_init(8, 10000, noise_variance)
            BgGrad.noise_variance = noise_variance
        else:
            assert self.noise_variance == noise_variance # Changing noise_variance currently requires restarting the process

    def to_cffi(self, img):
        # FIXME: assert continious
        return ffi.cast("uint8_t *", img.ctypes.data)
        
    def to_cffi16(self, img):
        return ffi.cast("uint16_t *", img.ctypes.data)

    def segment(self, img):
        assert img.dtype == np.uint8
        if self.bg is None:
            self.bg = np.zeros(img.shape, dtype=img.dtype)
            self.noise = np.zeros(img.shape, dtype=np.uint16)
            assert len(self.noise.shape) == 2
            self.noise[:,:] = int(5 * 256 / 1.48260221851)
            self.prev = None
            self.intensity = np.zeros(img.shape, dtype=img.dtype)
        if self.fcnt < self.maxfcnt:
            self.fcnt += 1
        grad = np.zeros(img.shape, dtype=img.dtype)
        step = 128.0 / self.fcnt
        istep = int(step)
        jump = 1
        if istep == 0:
            istep = 1
            jump = int(1.0 / step + 0.5)
            
        noise_step = 10 * 256 // self.fcnt

        if self.deshake:
            assert self.noise_variance is not ESTIMATE
            assert not self.exposure_check
            if self.blocksize == 1:
                lib.bggrad_deshake(self.to_cffi(img), img.shape[1], img.shape[0],
                                   self.to_cffi(grad), self.to_cffi(self.bg),
                                   istep, jump, self.deshake)
            else:
                lib.bggrad_block_deshake(self.to_cffi(img), img.shape[1], img.shape[0],
                                         self.to_cffi(grad), self.to_cffi(self.bg),
                                         istep, jump, self.deshake, self.blocksize)
        elif self.noise_variance is ESTIMATE:
            assert self.blocksize == 1
            prev = ffi.NULL if self.prev is None else self.to_cffi(self.prev)
            if self.exposure_check and self.fcnt >= self.maxfcnt:
                lib.bggrad_noise_black(self.to_cffi(img), img.shape[1], img.shape[0],
                                       self.to_cffi(grad), self.to_cffi(self.bg),
                                       istep, jump, prev, self.to_cffi16(self.noise),
                                       self.noise_scale, noise_step, self.to_cffi(self.intensity))
            else:
                lib.bggrad_noise(self.to_cffi(img), img.shape[1], img.shape[0],
                                 self.to_cffi(grad), self.to_cffi(self.bg),
                                 istep, jump, prev, self.to_cffi16(self.noise),
                                 self.noise_scale, noise_step)
            self.prev = img
        else:
            assert not self.exposure_check
            assert self.blocksize == 1
            lib.bggrad(self.to_cffi(img), img.shape[1], img.shape[0],
                       self.to_cffi(grad), self.to_cffi(self.bg),
                       istep, jump)
        return grad


if __name__ == '__main__':
    import sys
    assert sys.argv[1] is not None

    bg = BgGrad(blocksize=8, noise_variance=300, maxfcnt=5000) #blocksize=4, noise_variance=6**2, maxfcnt=100000)
    #bg = BgGrad(deshake=1, blocksize=16)
    
    import cv2
    vid_path = sys.argv[1]
    try:
        vid_path_int = int(vid_path)
        vid_path = vid_path_int
    except ValueError:
        pass
    vid = cv2.VideoCapture(vid_path)
    if not vid.isOpened():
        raise IOError(("Couldn't open video file or webcam. If you're "
        "trying to open a webcam, make sure your argument is an integer!"))
    
    accum_time = 0
    curr_fps = 0
    fps = "FPS: ??"
    prev_time = timer()
        
    while True:
        retval, image = vid.read()
        if not retval:
            print("Done!")
            sys.exit()

        # Convert to grayscale        
        grayscale = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        
        # Scale down a bit prior to running fgbg-segmentation
        grayscale = cv2.resize(grayscale, (1280,720))
        fg = bg.segment(grayscale)
                
        # Scale down again for visualization        
        to_show = np.vstack((fg, grayscale))
        to_show = cv2.resize(to_show, None, 0, 0.75, 0.75)
        
        curr_time = timer()
        exec_time = curr_time - prev_time
        prev_time = curr_time
        accum_time = accum_time + exec_time
        curr_fps = curr_fps + 1
        if accum_time > 1:
            accum_time = accum_time - 1
            fps = "FPS: " + str(curr_fps)
            curr_fps = 0
            
        cv2.rectangle(to_show, (0,0), (65, 22), (255,255,255), -1)
        cv2.putText(to_show, fps, (0,13), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0,0,0), 1)
                
        
        cv2.imshow("FGBG", to_show)
        cv2.waitKey(1)

