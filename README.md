# Three-Body-Problem-Gravitational-System
A python code to calculate planetary orbits in a three-body gravitational system. The code can demonstrate how one planet affects the orbit of another planet. As an example, Earth, Jupiter, Sun system is analyzed.

## Math and Theory
[Theory file](https://github.com/zaman13/Restricted-Three-Body-Problem-Gravitational-System/blob/master/Theory.ipynb)

## Requirements
The code uses matplotlib to create the orbital motion animation. The package ffmpeg is required to run the animation. It can be installed using the Anaconda terminal:
```python
conda install -c menpo ffmpeg
```
The video is embedded in the Jupyter notebook using html5. The vidoe can also exported as mp4. 



## Earth-Jupiter-Sun system
The effect of Jupiter on Earth's orbit can be visualized using the Earth-Jupiter-Sun system notebook. The mass of Jupiter was taken to be 500 times its original mass. The effect on Earth's orbit due to the presence of this "Super Jupiter" can be easily noticed from the orbital animation.


## Sample output
Earth-Jupiter-Sun system. 

<img src="https://github.com/zaman13/Restricted-Three-Body-Problem-Gravitational-System/blob/master/sample_output_1.png" alt="alt text" width="300">



## References
1. Earth fact sheet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
2. Jupiter fact sheet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
3. https://www.phas.ubc.ca/~berciu/TEACHING/PHYS349/karla.pdf
