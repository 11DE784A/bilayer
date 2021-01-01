# upscaled plot
upscale = 2.4
fntsm = Plots.font(pointsize=round(6.0*upscale))
fntlg = Plots.font(pointsize=round(8.0*upscale))
default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm, 
        legendtitlefont=fntlg)
default(size=(400*upscale,300*upscale)) #Plot canvas size
default(linewidth=upscale)
