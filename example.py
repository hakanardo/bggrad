from vi3o import Video, view, flipp

from bggrad import BgGrad
from bggrad.maxflow import MRFSegment

bg = BgGrad()
mrf = MRFSegment(0.05)

while True:
    for img in Video("test.mjpg", grey=True):
        fg = bg.segment(img)
        seg = mrf.segment(fg)

        flipp()
        view(img)
        view(fg)
        view(seg)