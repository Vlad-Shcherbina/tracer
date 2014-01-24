# Transform text from stdin in following format to bmp image
# ----------
# <filename>
# <width> <height>
# <r> <g> <b> ... <r> <g> <b>
# ...                    ...
# <r> <g> <b> ... <r> <g> <b>
# ----------

import PIL, Image, math

filename=raw_input()

s=raw_input()
s=s.split(" ")
width=int(s[0])
height=int(s[1])

im=Image.new("RGB",(width,height))

for j in range(height):
	s=raw_input()
	s=s.split(' ')
	for i in range(width):
		im.putpixel((i,j) , (int(s[i*3]),int(s[i*3+1]),int(s[i*3+2])))

#im.show()
im.save(filename,"BMP")
