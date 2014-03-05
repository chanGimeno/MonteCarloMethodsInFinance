#/usr/bin/env python

from pylab import *

##########################################################################################
def demo_eulerIntegration_bankAccount():
    """
    Demo_eulerIntegration_bankAccount: demo Euler integration scheme for exponential
    """
    ## parameters

    # intial conditions
    t0 = 0.0   # initial date
    M0 = 100.0 # initial deposit

    # interest rate
    r  = 0.05  

    # parameters for integration grid
    T  = 50.0  # length of time interval (in years)
    N  = 10     # number of integration steps

    ## Define integration grid
    deltaT = T/N                   # length of integration steps
    t      = linspace(t0,t0+T,N+1) # vector of monitoring times

    ## Euler integration method
    M_euler = zeros(N+1)
    M_euler[0] = M0  # initial condition

    for n  in range(N):
        M_euler[n+1] = M_euler[n]*(1 + r *deltaT)

    ## Exact solution
    M = M0*exp(r*(t-t0))

    ## Plot the results
    figure(1); clf()
    plot(t,M,'k-')
    plot(t,M_euler,'bo-')
    xlabel('t')
    ylabel('f(t)')
    legend(('M(t)',r'M$_{euler}$(t)'))

##########################################################################################
def eulerIntegration(t0,f0,a,T,N):
    """
    # eulerIntegration: Integration of an ODE with the Euler scheme
    
    # SYNTAX:
            [t,f] = eulerIntegration(t0,f0,a,T,N)
    
    # INPUT:
            t0 : Initial time
            f0 : Initial value of the function f(t0) = f0; 
            a  : Handle of the function a(t,f), which gives the value
                 of the derivative 
            T  : Length of integration interval [t0, t0+T]
            N  : Number of time steps
    
    # OUTPUT:
            t  : Times at which the trajectory is monitored
            f  : Values of the trajectory that starts from f(t0) = f0
                  
    
    # EXAMPLE:   
            # dsigma(t) = - alpha (sigma(t) -sigma_infty) dt   # reversion to the mean
            sigma_0 = 0.5; sigma_infty = 0.2; alpha = 1/2; 
            a = @(t,sigma)(-alpha*(sigma -sigma_infty))
            t0 = 0; T = 20;          
            N = 100;   
            [t,sigma] = eulerIntegration(t0,sigma_0,a,T,N);
            figure(1); plot(t,sigma);
            xlabel('t'); ylabel('f(t)'); legend('\sigma(t)');        
    """
    deltaT = 1.*T/N              # size of the integration step

    t = linspace(t0,t0+T,N+1) # monitoring times

    f = zeros(N+1)          # initialize trajectory

    # Euler integration method
    f[0] = f0                 # initial condition
    for n in range(N):
        f[n+1] = f[n] + a(t[n],f[n])*deltaT 

    return t,f

##########################################################################################
def demo_eulerIntegration_exponential():
    """
    Demo_eulerIntegration_exponential: Integration of ODE for the exponential  
    """

    # Exponential 
    # $$ dM(t) = r M(t) dt $$
    
    ## Parameters for the integration grid

    # initial conditions
    t0 = 0.  # initial time 
    M0 = 100. # initial value of the trajectory

    # length of integration interval
    T  = 50.

    ## Function that gives the value of the derivative 
    #
    # $$ a(t,M(t)) = r M(t) $$

    r = 0.05        # interest rate
    a = lambda t,M: r * M


    ## Integrate the ODE (coarse integration grid)    
    N = 10  # coarse grid
    [t, M_euler] = eulerIntegration(t0,M0,a,T,N)


    # # Plot the results

    # exact solution
    nPlot = 1000 
    tPlot = linspace(t0,t0+T,nPlot)   
    M     = M0 * exp( r * (tPlot-t0) )    # exact solution

    # compare the exact solution and the approximation by the Euler method 
    figure(1); clf()
    plot(tPlot,M,'k-')
    plot(t,M_euler,'bo-')
    xlabel('t'); ylabel('f(t)'); legend(('M(t)',r'M$_{euler}$(t)'))


    ## Integrate the ODE using a finer integration grid

    N = 50  # finer integration grid
    [t,M_euler] = eulerIntegration(t0,M0,a,T,N)

    ## Plot the results

    # compare the exact solution and the approximation given by the Euler method 
    figure(2); clf()
    plot(tPlot,M,'k-')
    plot(t,M_euler,'bo-')
    xlabel('t'); ylabel('f(t)'); legend(('M(t)',r'M$_{euler}$(t)'))

##########################################################################################
def demo_eulerIntegration_reversion2Mean():
    """
    # Demo_eulerIntegration_reversion2Mean: Integration of ODE with reversion to the mean  
    
    Reversion to the mean 
    
    $$ d\sigma(t) = - \alpha (\sigma(t) -\sigma_{\infty}) dt $$
    """

    ## Parameters for the integration grid

    T = 10.   # length of integration interval
    N = 100  # number of integration steps

    ## Function that gives the value of the derivative 
    #
    # $$ a(t,\sigma(t)) = - \alpha (\sigma(t) -\sigma_{\infty}) $$
    
    sigma_infty = 0.1 # reversion level
    alpha       = 1./2 # reversion rate

    a = lambda t,sigma: -alpha*(sigma - sigma_infty) 

    ## Integrate the ODE

    # initial conditions
    t0      = 0.   
    sigma_0 = 0.5 

    [t,sigma] = eulerIntegration(t0,sigma_0,a,T,N)
    
    # # Plot the results

    figure(1); clf()
    plotResults(t0,T,sigma_infty,t,sigma)

    ## Integrate the ODE with different initial conditions

    # initial conditions
    t0      = 0.  
    sigma_0 = 0.01 

    [t,sigma] = eulerIntegration(t0,sigma_0,a,T,N)

    # # Plot the results

    figure(2); clf()
    plotResults(t0,T,sigma_infty,t,sigma)

##########################################################################################
## Auxiliary function
def plotResults(t0,T,sigma_infty,t,sigma):
    """
    plot the results of demo_eulerIntegration_reversion2Mean
    """
    plot(t,sigma)
    reversionLevel = sigma_infty*ones(size(t))
    plot(t,reversionLevel,'k:')
    xlabel('t'); ylabel('f(t)') 
    legend((r'$\sigma$(t)',r'$\sigma_{\infty}$'),loc='best')
    xlim(t0, t0+T)
    ylim(0, max(sigma)*1.2)

##########################################################################################
def eulerIntegrationVariableGrid(f0,a,t):
    """
      eulerIntegrationVariableGrid: Integration of an ODE with the Euler scheme
     
      SYNTAX: 
             f = eulerIntegrationVariableGrid(f0,a,t)
     
      INPUT:
             f0   : Initial value of the function
             a    : Handle of the function a(t,f), which gives the value
                    of the derivative 
             t    : Times at which the trajectory is monitored
     
      OUTPUT:
             f    : Values of the trajectory that starts from 
                    f(1) = f0 at t(1) = t0
     
        EXAMPLE:   
             #  dM_t = r M_t dt   
             t = [1 3 7 8 10 12 15 20];
             M0 = 100; r = 0.05; a = @(t,M)(r*M);
             M_euler = eulerIntegrationVariableGrid(M0,a,t);
             nPlot = 1000; tPlot = linspace(t(1),t(end),nPlot); M = M0*exp(r*(tPlot-t(1)));
             figure(1); plot(tPlot,M,'k',t,M_euler,'b',t,M_euler,'bo','linewidth',2);
             xlabel('t'); ylabel('f(t)'); legend('M(t)','M_{euler}(t)',0)
    """

    deltaT = t[1:]-t[:-1]       # length of the different time intervals

    f = zeros(len(t)) # reserve memory for the trajectory

    f[0] = f0              # initial value of trajectory
    for n in range(len(deltaT)):
        f[n+1] = f[n] + a(t[n],f[n])*deltaT[n]

    return f

##########################################################################################
def demo_eulerIntegration_variableGrid():
    """
    demo_eulerIntegration_variableGrid: Integration of ODE for the exponential  
    """
    # Exponential 
    # $$ dM(t) = r M(t) dt $$

    # # Parameters for the integration grid

    # initial conditions
    M0 = 100. # initial value of the trajectory

    # Monitoring times for the integration
    t = array([1, 3, 7, 8, 10, 12, 15, 20])
         
    # # Function that gives the value of the derivative 
    #
    # $$ a(t,M(t)) = r M(t) $$

    r = 0.05        # interest rate
    a = lambda t,M: r * M

    # # Integrate the ODE (variable grid)

    M_euler = eulerIntegrationVariableGrid(M0,a,t)

    # # Plot the results

    # exact solution
    nPlot = 1000 
    tPlot = linspace(t[0],t[-1],nPlot)   
    M     = M0 * exp( r * (tPlot-t[0]) )    # exact solution

    # compare the exact solution and the approximation by the Euler method 
    figure(1); clf()
    plot(tPlot,M,'k-')
    plot(t,M_euler,'bo-')
    xlabel('t'); ylabel('f(t)'); legend(('M(t)',r'M$_{euler}$(t)'),loc='best')

##########################################################################################
def eulerIntegrationVanDerPol(t0,x0,p0,mu,T,N):
    """
    # eulerIntegrationVanDerPol: Van Der Pol oscillator solved via Euler integration
    
    # SYNTAX: 
            [t,x,p] = eulerIntegrationVanDerPol(t0,x0,p0,mu,T,N)
    
    # INPUT:
            t0 : Initial time
            x0 : Initial position 
            p0 : Initial momentum 
            mu : damping strength
            T  : Length of integration interval [t0, t0+T]
            N  : Number of time steps
    
    # OUTPUT:
            t  : Times at which the trajectory is monitored
                 t(n) = t0 + n Delta T
            x  : Values of the position along the trajectory
            p  : Values of the momentum along the trajectory
                 
    # EXAMPLE 1:   
            t0 = 0; x0 = 0; p0 = 1;
            mu = 0.01;
            T  = 20; N = 1e5;   
            [t,x,p] = eulerIntegrationVanDerPol(t0,x0,p0,mu,T,N);
            figure(1); plot(t,x,t,p);
            xlabel('t');legend('x(t)','p(t)');        
              
    # EXAMPLE 2:   
            t0 = 0; x0 = 0; p0 = 1;
            mu = 2;
            T = 20; N = 1e5;   
            [t,x,p] = eulerIntegrationVanDerPol(t0,x0,p0,mu,T,N);
            figure(1); plot(t,x,t,p);
            xlabel('t');legend('x(t)','p(t)');        
    """
    deltaT = 1.*T/N              # size of integration step
    
    t = linspace(t0,t0+T,N+1) # initialize monitoring times
    
    x = zeros(N+1)          # initialize x
    p = zeros(N+1)          # initialize p
    
    # Euler integration 
    x[0] = x0                 # initial conditions
    p[0] = p0                 
    for n in range(N):
        x[n+1] = x[n] + p[n]*deltaT 
        p[n+1] = p[n] + (mu*(1-x[n]*x[n])*p[n]-x[n])*deltaT 

    return t, x, p

##########################################################################################
def demo_vanDerPol():
    """
    demo_vanDerPol: numerical integration (Euler) of the Van der Pol oscillator
    """
    # damping parameter
    mu = 2

    # initial conditions
    t0 = 0. 
    x0 = 1e-5 
    p0 = 1e-5

    # integration parameters
    T = 100. 
    N = int(1e5)

    [t,x,p] = eulerIntegrationVanDerPol(t0,x0,p0,mu,T,N)

    # Plot  p(t) as a function of x(t)
    f1 = figure(1); clf()
    h1, = plot(x0,p0,'k',linewidth=1.5)
    h2, = plot(x0,p0,'ro',linewidth=2,markersize=6)
    MIN_X = -2.2;  MAX_X = 2.2
    MIN_Y = -4;    MAX_Y = 4
    xlim(MIN_X, MAX_X)
    ylim(MIN_Y, MAX_Y)
    xlabel('x(t)')
    ylabel('p(t)')
    draw()

    # Plot time series of x(t) and p(t)
    f2 = figure(2); clf()
    subplot(2,1,1)
    h3, = plot(t0,x0)
    xlim(t0, t0+T)
    ylim(MIN_X, MAX_X)
    xlabel('t'); ylabel('x(t)')
    subplot(2,1,2)
    h4, = plot(t0,p0,'m')
    xlim(t0, t0+T)
    ylim(MIN_Y, MAX_Y)
    xlabel('t'); ylabel('p(t)')
    draw()

    pause(3)

    # Update plots
    frameLength = 500
    for i in range(len(t)):
        if(mod(i+1,frameLength) == 0):
            h1.set_data(x[:i], p[:i])
            h2.set_data(x[i], p[i])
            f1.canvas.draw()
            h3.set_data(t[:i], x[:i])
            h4.set_data(t[:i], p[:i])
            # pause(0.1)
            f2.canvas.draw()

##########################################################################################
if __name__=="__main__":
    
    # demo bank account
    close("all")
    demo_eulerIntegration_bankAccount()
    pause(5)

    # demo integration exponential
    close("all")
    demo_eulerIntegration_exponential()
    pause(5)
         
    # demo_eulerIntegration_reversion2Mean()
    close("all")
    demo_eulerIntegration_reversion2Mean()
    pause(5)
         
    # demo_eulerIntegration_variableGrid
    close("all")
    demo_eulerIntegration_variableGrid()
    pause(5)

    # demo_vanDerPol
    close('all')
    demo_vanDerPol()
