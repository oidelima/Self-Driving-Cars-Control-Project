function [df_dx, df_dy] = brianNumGradient(f,x,y,step_size)
    df_dx = (f(x+step_size,y) - f(x,y)) / step_size;
    df_dx = df_dx + ((f(x,y) - f(x-step_size,y)) / step_size);
    df_dx = df_dx / 2;
    
    df_dy = (f(x,y+step_size) - f(x,y)) / step_size;
    df_dy = df_dy + ((f(x,y) - f(x,y-step_size)) / step_size);
    df_dy = df_dy / 2;
end

