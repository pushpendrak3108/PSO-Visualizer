function fitness = GetFitness( position )
% input coordinates x and y output fitness value
% function value z = fitness value

x=position(1);
y=position(2);
	  
fitness = FuncCalculate(x, y);

end

