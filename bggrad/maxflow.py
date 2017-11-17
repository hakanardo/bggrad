from array import array
from bggrad._maxflow import ffi, lib
from math import log
import numpy as np

def void_p(data):
    if isinstance(data, np.ndarray):
        p, _ = data.__array_interface__['data']
    else:
        p, _ = data.buffer_info()
    return ffi.cast('void *', p)

class Graph(object):
    def __init__(self, nodes):
        self.nodes = nodes
        self.graph = lib.new_graph(nodes, 2*nodes)
        self.reuse = False
        self.old_source = array('i', [0]) * nodes
        self.old_sink = array('i', [0]) * nodes
        
    def set_capacity(self, source, sink):
        assert len(source) == len(sink)
        dsource = array('i', (c-o for c,o in zip(source, self.old_source)))
        dsink = array('i', (c-o for c,o in zip(sink, self.old_sink)))
        lib.add_capacity(self.graph, len(source), void_p(dsource), void_p(dsink))
        self.old_source[:] = array('i', source)
        self.old_sink[:] = array('i', sink)

    def set_capacity_from_image(self, prior, image, scale):
        assert image.dtype == np.uint8 and len(image.shape) == 2
        assert len(image.flat) == self.nodes
        lib.set_capacity_from_image(self.graph, self.nodes, void_p(self.old_source), void_p(self.old_sink), 
                                    void_p(image), prior, scale)
        
    def add_edge(self, head, tail, cap):
        self.reuse = False
        lib.add_edge(self.graph, head, tail, cap)
        
    def add_edges(self, heads, tails, caps):
        heads, tails, caps = array('i', heads), array('i', tails), array('i', caps)
        assert len(heads) == len(tails) == len(caps)
        lib.add_edges(self.graph, len(heads), void_p(heads), void_p(tails), void_p(caps))
        
    def maxflow(self):
        flow = lib.maxflow(self.graph, self.reuse)
        #self.reuse = True
        return flow
        
    def segments(self):
        a = array('i', [-1]) * self.nodes
        lib.get_segments(self.graph, self.nodes, void_p(a))
        return a

        
class MRFSegment(object):
    def __init__(self, pdiff, prior=0.5):
        self.pdiff = pdiff
        self.prior = prior
        self.graph = None
        self.scale = 1000
        self.eps = 0.0001
    
    def segment(self, img):
        if self.graph is None:
            height, width = img.shape
            self.graph = Graph(width * height)
            heads, tails = [], []
            for i, (y, x) in enumerate(np.ndindex(*img.shape)):
                if x > 0:
                    heads.append(i-1)
                    tails.append(i)
                if y > 0:
                    heads.append(i-width)
                    tails.append(i)
                if 0 < x < width-1 and y > 0:
                    heads.append(i-width-1)
                    tails.append(i)
                    heads.append(i-width+1)
                    tails.append(i)
                    
            cap = int((-log(self.pdiff + self.eps) - -log(1.0 - self.pdiff + self.eps)) * self.scale)
            self.graph.add_edges(heads, tails, [cap] * len(heads))
        
        self.graph.set_capacity_from_image(self.prior, img, self.scale)

        self.graph.maxflow()
        seg = self.graph.segments()
        return (255 * np.array(seg).reshape(img.shape)) #.astype(np.uint8)
        # res = img.new(typecode='B')
        # for i in range(len(img.flat)):
        #     res.flat[i] = seg[i] * 255
        # return res

            
            
        
