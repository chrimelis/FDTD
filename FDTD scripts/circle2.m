function circle2(x,y,r)
    % used for plotting a MAGENTA circle
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(xunit, yunit, 'm','LineWidth',2);
end