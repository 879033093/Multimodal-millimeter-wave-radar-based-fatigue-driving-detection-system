function [retVal1, retVal2, retVal3] = PlotsCentroid(rx_2dcfar)

rx_2dcfar_temp = padarray(rx_2dcfar, [1,1], 0);
[r, c] = size(rx_2dcfar);
num_target = 0;
add = zeros(1000, 2);

for i = 1:r
    for j = 1:c
        if rx_2dcfar_temp(i+1, j+1) > 0
            a = rx_2dcfar_temp(i:i+2, j:j+2);
            b = max(max(a));
            [x, y] = find(a == max(max(a)));
            temp = zeros(3, 3);
            rx_2dcfar_temp(i:i+2, j:j+2) = temp;
            temp(x, y) = b;
            rx_2dcfar_temp(i:i+2, j:j+2) = temp;
        end
    end
end

rx_2dcfar_plots = rx_2dcfar_temp(2:r+1, 2:c+1);

for i = 1:r
    for j = 1:c
        if rx_2dcfar_plots(i, j) > 0
            num_target = num_target + 1;
            add(num_target, 1) = i;
            add(num_target, 2) = j;
        end
    end
end

add = add(1:num_target, :);

retVal1 = rx_2dcfar_plots;
retVal2 = add;
retVal3 = num_target;
