#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 15:32:58 2019

@author: jloos
"""

def plot(t):
    
    Plot_hyd
    
    
    
    
    
    
    
    
    fig.canvas.draw()       # draw the canvas, cache the renderer
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return image


kwargs_write = {'fps':4.0, 'quantizer':'nq'}
imageio.mimsave('./powers.gif', [plot(t) for t in range(200)], fps=4)