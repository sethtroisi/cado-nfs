

# Maybe it's just me, but frankly the native sage plot interface sucks.
# This is a wrapper around gnuplot, requiring the gnuplotpy sage package.
# (currently gnuplotpy-1.7.p3.spkg, works ok for me).

# This program does spit out a lot of mess in /tmp/

# In fact, we don't even use Gnuplot.PlotItems ; it's much easier the
# straight way, without pipes or fancy tricks.
# import Gnuplot.PlotItems


#
#def myplot_list(foo, fcn, x0, x1, np):
#    # It's just so incredibly slow...
#    t=(x1-x0)/np
#    x=x0
#    a=[ [] for i in range(1 + len(fcn))]
#    for i in [0..np]:
#        a[0].append(float(x))
#        for j in range(len(fcn)):
#            a[1+j].append(float(fcn[j](x)))
#        x+=t

def myplot_list(fcn0, x0, x1, np):
    t=(x1-x0)/np
    fcn=[ z if type(z) == tuple else (z,z.func_name) for z in fcn0 ]
    h=reduce(lambda x,y: hash((x,y)), fcn, hash((x0,x1)))
    filename="/tmp/gnuplot.%u" % abs(h)
    f=open(filename,"w")
    x=x0
    for i in [0..np]:
        s=reduce(operator.add, [" %f"  % float(g[0](x)) for g in fcn])
        f.write("%f%s\n" % (x,s))
        x+=t
    f.close()
    time.sleep(0.1)
    pattern="'%s' using 1:%d with lines title '%s'"
    command="plot "
    for i in range(len(fcn)):
        if i > 0: command+=", "
        command += pattern % (filename, i+2, fcn[i][1])
    gnuplot.gnuplot().clear()
    gnuplot.gnuplot().reset()
    gnuplot.gnuplot().gnuplot("set terminal wxt\n")
    gnuplot.gnuplot().gnuplot(command)
    # pl=Gnuplot.PlotItems.File(filename, **{'with': 'lines'})
    # foo.plot(pl)
    # We no longer unlink the file. I know, that's dirty.

def myplot(fcn, x0, x1, np=200):
    """plots the given function(s) between x0 and x1. The fcn argument
    may be either a single function, or a list of such. In turn,
    functions may be lambda terms, or tuples (lambda, name), so that the
    given name is printed onthe plot"""
    if type(fcn) == list:
        myplot_list(fcn, x0, x1,np)
    else:
        myplot_list([fcn],x0,x1,np)


