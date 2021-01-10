import os
from PIL import Image
import numpy as np
import requests
from io import BytesIO

def map_url(bodyName, mapType):
    mapDirURL = "https://raw.githubusercontent.com/theastrogoth/KSP-Trajectory-Illustrator/master/assets/images/"

    return mapDirURL + str(bodyName) + str(mapType) + '.png'

def combine_tiles(path=None):
    
    if not path is None:
        os.chdir(path)
    
    for path, directories, files in os.walk(path):
        if 'Tile0000.png' in files and 'Tile0001.png' in files:
    
            images = [Image.open(x) for x in [os.path.join(path,'Tile0000.png'), os.path.join(path,'Tile0001.png')]]
            widths, heights = zip(*(i.size for i in images))
            
            total_width = sum(widths)
            max_height = max(heights)
            
            new_im = Image.new('RGB', (total_width, max_height))
            
            x_offset = 0
            for im in images:
              new_im.paste(im, (x_offset,0))
              x_offset += im.size[0]
            
            new_im.save(os.path.join(path,'combined.png'))

def make_small_image(path=None):
    
    if not path is None:
        os.chdir(path)
    
    for path, directories, files in os.walk(path):
        for file in files:
            if file[-4:] == '.png':
                image = Image.open(os.path.join(path, file))
                width, height = image.size
                
                if width==2048 and height==1024:
                    image = image.convert('RGB')
                    newImage = image.resize((512, 512))
                    pix = image_colormap(list(newImage.getdata()))[2]
                    newImage.putdata(pix)
                    newImage.save(os.path.join(path,file[:-4]+'Small.png'))

def round_colors(pix, roundVal):
    
    newPix = []
    for pick in pix:
        rgb = list(pick)
        for jj, cl in enumerate(pick):
            # rgb[jj] = cl - cl%roundVal + roundVal/2
            rgb[jj] = round(roundVal * round(cl/roundVal))
        newPix.append(tuple(rgb))
    
    return newPix

def get_pixel_values(imagePath, url=True):
    if url:
        response = requests.get(imagePath)
        image = Image.open(BytesIO(response.content))
    else:
        image = Image.open(imagePath)
    image = image.convert('RGB')
    pix = list(image.getdata())
    width, height = (image.size)
    return pix, width, height

def pixel_to_grayscale_value(r=None,g=None,b=None, tup=None):
    
    if not tup is None:
        r = tup[0]
        g = tup[1]
        b = tup[2]
    
    return (r+b+g)/3

def image_colormap(pix, roundNum=85, rounded=False):
    
    if not rounded:
        numUnique = 999
        # plotly seems to fail if the colormap has >~100 intervals
        while numUnique > 100:
            newPix = round_colors(pix, 255/roundNum)
            
            uniqueColors = np.unique(newPix, axis=0)
            numUnique = len(uniqueColors)
            if roundNum <= 5:
                roundNum = roundNum-1
            else:
                roundNum = roundNum - round(roundNum/10)
    else:
        newPix = pix
        uniqueColors = np.unique(pix, axis=0)
        numUnique = len(uniqueColors)
    
    breaks = np.linspace(0,1,numUnique+1)
    interval = breaks[1]
    mapVal = 0
    
    colorMap = [[0, 'rgb'+str(tuple(uniqueColors[0]))]]
    mapDict = {}
    mapDict[tuple(uniqueColors[0])] = interval/2
    for ii in range(numUnique):
        colorMap.append([mapVal, 'rgb'+str(tuple(uniqueColors[ii]))])
        mapDict[tuple(uniqueColors[ii])] = mapVal+interval/2
        mapVal = mapVal + interval
        if ii == numUnique-1:
            mapVal = 1
        colorMap.append([mapVal, 'rgb'+str(tuple(uniqueColors[ii]))])
        
    return colorMap, mapDict, newPix
