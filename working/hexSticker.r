### mcmcHammer hex sticker
###
### source('E:/Ecology/Drive/R/enmSdmX/working/hexSticker.r')

library(hexSticker)
library(magick)
library(sysfonts)

img <- image_read('E:/Ecology/Drive/R/mcmcHammer/working/enmSdmIcon.png')

sticker(
	subplot = img,
	package='enmSdmX',
	p_size=18,
	p_color='gray90',
	# p_y=1.4,
	p_y=0.6,
	s_x=1,
	# s_y=0.75,
	s_y=1.25,
	s_width=0.97,
	s_height=0.97,
	h_fill = 'black',
	white_around_sticker = FALSE,
	filename='E:/Ecology/Drive/R/mcmcHammer/working/enmSdm.png'
)


		