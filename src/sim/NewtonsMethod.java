/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sim;

import static java.lang.Math.abs;

/**
 *
 * @author Logan
 */
public class NewtonsMethod {
    public double target;
    public double x1, x2;
    public double y1, y2;
    public double epsilon;
    public double sgn;
    public boolean done = false;
    public double finalX;
    public double finalY;
    public boolean error=false;
    public RealFunc F;
    public int stepcount=2;
    public int stct=0;
    public double[] xSteps=new double[100];
    public double[] x1Steps=new double[100];
    public double[] x2Steps=new double[100];
    public double[] y1Steps=new double[100];
    public double[] y2Steps=new double[100];
    public double[] ySteps=new double[100];

    public NewtonsMethod(double target, double x1, double x2, double y1, double y2, double epsilon){
        this.target=target;
        this.x1=x1;
        this.x2=x2;
        this.y1=y1;
        this.y2=y2;
        this.epsilon=epsilon;
        xSteps[0]=x1;
        xSteps[1]=x2;
        ySteps[0]=y1;
        ySteps[1]=y2;
               
        
        if (abs(y1 - target) < epsilon){
            //done
            done = true;
            finalX = x1;
            finalY = y1;
        } else if (abs(y2 - target) < epsilon) {
            //done            
            done = true;
            finalX = x2;
            finalY = y2;
        }
    }
    
    public double nextStepX(){
        
        x1Steps[stct]=x1;
        x2Steps[stct]=x2;
        y1Steps[stct]=y1;
        y2Steps[stct]=y2; 
        stct++;
        return x1+(target - y1) * (x2-x1) / (y2 - y1);
    }
    
    public void nextStep(double x, double y){
        xSteps[stepcount]=x;
        ySteps[stepcount]=y;
        stepcount++;
        if (abs(y - target) < epsilon){
            //done
            done = true;
            finalX = x;
            finalY = y;
        } else {
            x1=x2;
            y1=y2;
            x2=x;
            y2=y;            
        }
    }
    
    public double solveForX(){
        while (!done){
            double x = nextStepX();
            double y = F.eval(x);
            nextStep(x,y);
        }
        
        return finalX;
    }
}
