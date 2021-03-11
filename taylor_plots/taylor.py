import numpy
import matplotlib.artist as artist
import matplotlib.patches as patches
import matplotlib.lines as lines
import matplotlib.pylab as pylab
from matplotlib.axes import Axes

def regrid(data1,lat1,lon1,lat2,lon2):

    # check inputs
    if len(lat1) < len(lat2):
        print 'Target grid is higher resolution than input grid.'
        raise
    if len(lon1) < len(lon2):
        print 'Target grid is higher resolution than input grid.'
        raise

    shape = (len(lat2),len(lon2))
    data2 = numpy.ma.zeros(shape,dtype=data1.dtype,)
    if 'fill_value' in dir(data1):
        data2.fill_value = data1.fill_value

    # regrid data
    for idx2,lon in enumerate(lon2):
        for idy2,lat in enumerate(lat2):
            idx1 = numpy.abs(lon1 - lon).argmin()
            idy1 = numpy.abs(lat1 - lat).argmin()
            data2[...,idy2,idx2] = data1[...,idy1,idx1]
    return data2


class Taylor_statistics():

    cc = 0.0
    ratio = 0.0
    bias = 0.0
    rmsnorm = 0.0

    def __init__(
            self,
            test_data,test_lat,test_lon,
            cntl_data,cntl_lat,cntl_lon,
        ):

        # regrid data
        if len(test_lat) < len(cntl_lat):
            lat = test_lat
        else:
            lat = cntl_lat
        if len(test_lon) < len(cntl_lon):
            lon = test_lon
        else:
            lon = cntl_lon
        test_data = regrid(test_data,test_lat,test_lon,lat,lon)
        cntl_data = regrid(cntl_data,cntl_lat,cntl_lon,lat,lon)

        # calculate weights
        lat_weights = numpy.cos(lat[:]*numpy.pi/180.)
        lat_weights = lat_weights/lat_weights.sum()

        # make full lat,lon dimensioned weight array
        weight_array = numpy.ma.ones([len(lat),len(lon)])
        for idx in range(len(lon)):
            weight_array[:,idx] = lat_weights[:]

        # mask missing values
        test_data.mask = numpy.ma.mask_or(test_data.mask,cntl_data.mask)
        cntl_data.mask = numpy.ma.mask_or(test_data.mask,cntl_data.mask)
        weight_array.mask = test_data.mask

        # calculate statistics
        self.calculate(test_data,cntl_data,weight_array)

    def calculate(self,test,cntl,wgt):
        """Calculate Taylor statistics for making taylor diagrams."""

        # calculate sums and means
        sumwgt = numpy.ma.sum(wgt)
        meantest = numpy.ma.sum(wgt*test)/sumwgt
        meancntl = numpy.ma.sum(wgt*cntl)/sumwgt

        # calculate variances
        stdtest = (numpy.ma.sum(wgt*(test-meantest)**2.0)/sumwgt)**0.5
        stdcntl = (numpy.ma.sum(wgt*(cntl-meancntl)**2.0)/sumwgt)**0.5

        # calculate correlation coefficient
        ccnum = numpy.ma.sum(wgt*(test-meantest)*(cntl-meancntl))
        ccdem = sumwgt*stdtest*stdcntl
        self.cc = ccnum/ccdem

        # calculate variance ratio
        self.ratio = stdtest/stdcntl

        # calculate bias
        self.bias = (meantest - meancntl)/numpy.abs(meancntl)
        #self.bias = meantest - meancntl

        # calculate centered pattern RMS difference
        rmssum = numpy.ma.sum(wgt*((test-meantest)-(cntl-meancntl))**2.0)
        rmserr = (rmssum/sumwgt)**0.5
        self.rmsnorm = rmserr/stdcntl


class Taylor_diagram():
    """
    Custom axes class to facilitate making Taylor diagrams. Adds the taylor
    method to draw the diagram.
    """

    def __init__(
            self,ax,cc,ratio,bias,
            rms=None,casecolors=None,varlabels=None,
        ):

        # default casecolors
        if casecolors is None:
            casecolors = [
                    'blue','red','green','cyan','magenta','yellow',
                    '0.75','0.50','0.25','0.10',
                ]

        # copy axes instance
        self.ax = ax

        # plot size
        #self.xymax = 1.50
        self.xymax = numpy.max([1.50,numpy.max(ratio)+numpy.max(bias)/2.0])

        # draw axes
        self.draw_axes(ratio,cc)

        # add some reference arcs
        self.draw_ratio(1.0,linestyle='solid',color='black')
        self.draw_ratio(numpy.min(ratio),linestyle='dotted',color='black')
        self.draw_ratio(numpy.max(ratio),linestyle='dotted',color='black')

        # draw some reference lines
        self.draw_radii(numpy.min(cc))
        self.draw_radii(numpy.max(cc))

        # draw points
        self.draw_point(
                ratio,cc,bias,
                bubblecolors=casecolors,
                varlabels=varlabels,
            )

        # draw rms circles
        if rms is not None:
            self.draw_rms(min(rms))
            self.draw_rms(max(rms))

    def polar_transform(self,r,theta):
        x = r*numpy.cos(theta)
        y = r*numpy.sin(theta)
        return x,y

    def draw_radii(self,cc,color='black',linestyle='dotted'):
        theta = numpy.arccos(cc)
        x,y = self.polar_transform(self.xymax,theta)
        self.ax.plot([0,x],[0,y],':k',linewidth=0.5)

    def draw_ratio(self,r,color='black',linestyle='dotted'):
        arc = patches.Arc(
                [0,0],2*r,2*r,
                theta1=0.0,
                theta2=90.0,
                color=color,
                linestyle=linestyle,
                linewidth=0.5,
            )
        self.ax.add_patch(arc)

    def draw_axes(self,ratio,cc):

        # draw axes
        wedge = patches.Wedge(
                [0,0],self.xymax,
                theta1=0.0,
                theta2=90.0,
                color='black',
                fill=False
            )
        self.ax.add_patch(wedge)

        # correlation axis label
        self.ax.text(
                self.xymax*0.8,self.xymax*0.8,'Correlation',
                color='black',
                rotation='-45',
                horizontalalignment='center',
                verticalalignment='center'
            )

        # x-axis label
        self.ax.xaxis.set_label_text('Relative standard deviation')

        # x-axis tickmarks
        xmajorticks = numpy.arange(0.0,self.xymax+0.01,0.25)
        xmajorticklabels = []
        xminorticks = [numpy.min(ratio),numpy.max(ratio)]
        xminorticklabels = ["%.2f"%(numpy.min(ratio)),"%.2f"%(numpy.max(ratio))]

        # y-axis tickmarks
        ymajorticks = numpy.arange(0.0,self.xymax+0.01,0.25)
        ymajorticklabels = []
        for idy in range(len(ymajorticks)):
            if idy%2 == 0:
                ymajorticklabels.append('%.1f'%(ymajorticks[idy]))
            else:
                ymajorticklabels.append('')

        # correlation tickmarks
        ccmajorticks = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        ccminorticks = [0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99]
        ccmajorticklabellocs = [
                0.5,0.95,numpy.min(cc),numpy.max(cc)
            ]
        ccmajorticklabels = [
                "0.5","0.95","%.2f"%(numpy.min(cc)),"%.2f"%(numpy.max(cc))
            ]

        # set and draw tickmarks
        self.ax.set_aspect('equal')
        self.ax.xaxis.set_ticks_position('bottom')
        self.ax.yaxis.set_ticks_position('left')

        self.ax.xaxis.set_ticks(xmajorticks)
        self.ax.xaxis.set_ticklabels(xmajorticklabels)
        self.ax.xaxis.set_ticks(xminorticks,minor=True)
        self.ax.xaxis.set_ticklabels(xminorticklabels,minor=True)

        self.ax.yaxis.set_ticks(ymajorticks)
        self.ax.yaxis.set_ticklabels(ymajorticklabels)


        # draw correlation tickmarks
        majorticklength = self.xymax*0.02
        minorticklength = self.xymax*0.01
        for cctick in ccmajorticks:
            r = [self.xymax-majorticklength,self.xymax            ]
            theta = [numpy.arccos(cctick),numpy.arccos(cctick)]
            x,y = self.polar_transform(r,theta)
            self.ax.plot(x,y,'-k')
        for i,ccticklabel in enumerate(ccmajorticklabels):
            x,y = self.polar_transform(self.xymax,numpy.arccos(ccmajorticklabellocs[i]))
            self.ax.text(x,y,ccticklabel,color="black")
        for cctick in ccminorticks:
            r = [self.xymax-minorticklength,self.xymax            ]
            theta = [numpy.arccos(cctick),numpy.arccos(cctick)]
            x,y = self.polar_transform(r,theta)
            self.ax.plot(x,y,'-k')
 
        # re-draw xy ticks to be consistent with cc ticks (not elegant)
        for xx in xmajorticks:
            x = [xx,xx                         ]
            y = [0 ,majorticklength]
            self.ax.plot(x,y,'-k')
        for yy in ymajorticks:
            x = [0 ,majorticklength]
            y = [yy,yy                         ]
            self.ax.plot(x,y,'-k')


    def draw_point(self,ratio,cc,bias,bubblecolors=None,varlabels=None,labeloffset=0.025):
        
        if len(ratio.shape) == 2:
            nvars = len(ratio[:,0])
            ncases = len(ratio[0,:])
            for ivar in range(nvars):

                # transform coordinates
                r = ratio[ivar,:]
                theta = numpy.arccos(cc[ivar,:])
                size = numpy.abs(bias[ivar,:])/2.0
                x,y = self.polar_transform(r,theta)
                for icase in range(ncases):

                    # draw "bubbles" around points
                    if bubblecolors is None: 
                        bubblecolor='blue'
                    else:
                        bubblecolor=bubblecolors[icase]

                    circle = patches.Circle(
                            (x[icase],y[icase]),size[icase],
                            color=bubblecolor,
                            linewidth=0.0,
                            alpha=0.30,
                        )
                    self.ax.add_patch(circle)

                    # draw the actual points
                    pointcolor=bubblecolor
                    circle = patches.Circle(
                            (x[icase],y[icase]),0.01,
                            color=pointcolor,
                        )
                    self.ax.add_patch(circle)

                    # add labels to points if given
                    if varlabels is not None:
                        self.ax.text(
                                x[icase],y[icase]+labeloffset,
                                varlabels[ivar],
                                horizontalalignment='center',
                                verticalalignment='bottom',
                                color=bubblecolor,
                            )
        else:

            # transform coordinates
            r = ratio
            theta = numpy.arccos(cc)
            size = numpy.abs(bias)/2.0
            x,y = self.polar_transform(r,theta)
            for icase in range(len(x)):

                # draw "bubbles" around points
                if bubblecolors is None: 
                    bubblecolor='blue'
                else:
                    bubblecolor=bubblecolors[icase]

                circle = patches.Circle(
                        (x[icase],y[icase]),size[icase],
                        color=bubblecolor,
                        linewidth=0.0,
                        alpha=0.30,
                    )
                self.ax.add_patch(circle)

                # draw the actual points
                pointcolor=bubblecolor
                circle = patches.Circle(
                        (x[icase],y[icase]),0.01,
                        color=pointcolor,
                        #edgecolor="black",
                        #linewidth=0.5,
                    )
                self.ax.add_patch(circle)

                # add label to points if given
                if varlabels is not None:
                    self.ax.text(
                            x[icase],y[icase]+labeloffset,labels[icase],
                            horizontalalignment='center',
                            verticalalignment='bottom',
                            color=bubblecolor,
                        )

    def draw_rms(self,r,color='black',linestyle='dotted'):
        tmp = numpy.pi - numpy.arccos((r**2.0 + 1 - self.xymax**2.0)/2*r)
        theta1 = numpy.max([0.0,180.0/numpy.pi*tmp])
        theta2 = 180.0
        arc = patches.Arc(
                [1,0],2*r,2*r,
                theta1=theta1,
                theta2=theta2,
                color=color,
                linestyle=linestyle,
                linewidth=0.5,
            )
        self.ax.add_patch(arc)
