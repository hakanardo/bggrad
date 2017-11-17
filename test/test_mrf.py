from random import randrange
from vi3o import flipp, view, viewsc
from bggrad.maxflow import Graph, MRFSegment
import numpy as np


class TestMRF:
    def test_graph(self):
        g = Graph(2)
        g.set_capacity([2, 2], [0, 6])
        g.add_edge(0, 1, 1)
        flow = g.maxflow()
        assert flow == 3
        assert tuple(g.segments()) == (1, 0)
        assert g.maxflow() == 3
        assert tuple(g.segments()) == (1, 0)
        g.set_capacity([2, 2], [6, 0])
        flow = g.maxflow()
        assert flow == 3
        assert tuple(g.segments()) == (0, 1)

    def test_segment(self):
        segmenter = MRFSegment(0.05)

        img = np.zeros([240, 320], np.uint8)
        img[10:50, 30:60] = 200
        for i in range(100):
            x, y = randrange(img.shape[1]), randrange(img.shape[0])
            img[y, x] = 255 - img[y, x]
        fg = segmenter.segment(img)
        if False:
            flipp(pause=True)
            view(img)
            viewsc(fg)
        assert fg[42, 42] == 255
        assert fg[142, 142] == 0
