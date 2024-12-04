function [cooling_rate_Post] = POST_cooling_carbon(temperature)
%POST_CARBON Summary of this function goes here

cooling_rate_Post = 0.0 * temperature;
for i = 1: length(temperature)
   if( temperature(i) < 3.0 )
%        cooling_rate_Post(i) =  -1.901213e1;
       cooling_rate_Post(i) = ( 1.965300e3 +log10(temperature(i)/1.0e3)*( 4.572039e3 +log10(temperature(i)/1.0e3)*( 4.159590e3 +log10(temperature(i)/1.0e3)*( 1.871560e3 +log10(temperature(i)/1.0e3)*( 4.173889e2 +log10(temperature(i)/1.0e3)* 3.699382e1 )))));
   elseif( temperature(i) < 20.0 )
       cooling_rate_Post(i) = ( 1.965300e3 +log10(temperature(i)/1.0e3)*( 4.572039e3 +log10(temperature(i)/1.0e3)*( 4.159590e3 +log10(temperature(i)/1.0e3)*( 1.871560e3 +log10(temperature(i)/1.0e3)*( 4.173889e2 +log10(temperature(i)/1.0e3)* 3.699382e1 )))));
   elseif( temperature(i) < 200.0 )
       cooling_rate_Post(i) = ( 7.467599e1 +log10(temperature(i)/1.0e3)*( 4.549038e2 +log10(temperature(i)/1.0e3)*( 8.372937e2 +log10(temperature(i)/1.0e3)*( 7.402515e2 +log10(temperature(i)/1.0e3)*( 3.147607e2 +log10(temperature(i)/1.0e3)* 5.164578e1 )))));
   elseif( temperature(i) < 2000.0 )
       cooling_rate_Post(i) = (-2.120151e1 +log10(temperature(i)/1.0e3)*(-3.668933e-1+log10(temperature(i)/1.0e3)*( 7.295099e-1+log10(temperature(i)/1.0e3)*(-1.944827e-1+log10(temperature(i)/1.0e3)*(-1.263576e-1-log10(temperature(i)/1.0e3)* 1.491027e-1)))));
   end
end
cooling_rate_Post = (1.e-13/1.602)*(10.0).^(19.0+cooling_rate_Post);

end

