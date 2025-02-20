function number_array = get_coeffs(file_name,type,reaction)
    marker = 'BEGIN DATA HERE';
    begin_table = '\begin{small}\begin{verbatim}';
    % Open the file for reading
    fileID = fopen(file_name, 'r');

    if fileID == -1
        error('Unable to open the file.');
    end

    found_marker = false;
    found_type = false;
    found_reaction = false;
    found_begin_table = false;
    string = '';
    line_nr = 0;
    line_counter = 0;

    % Read the file line by line
    while ~feof(fileID)
        line_nr = line_nr+1;

        % Read the current line
        line = fgetl(fileID);

        % Check if the marker is found
        if ~found_marker
            if contains(line, marker)
                found_marker = true;
%                 fprintf('Marker found at line %d\n', line_nr);
            end
        else
            if  ~found_type
                if contains(line, append('\section{',type))
                    found_type = true;
%                     fprintf('Type found at line %d\n', line_nr);
                end
            else
                if  ~found_reaction
                    l = strip(line);
                    s = append('Reaction ', reaction, ' ');
                    if length(l)>=length(s)
                        if strcmp(l(1:length(s)), s)
                            found_reaction = true;
    %                         fprintf('Reaction found at line %d\n', line_nr);
                        end
                    end
                else
                    if ~found_begin_table
                        if contains(line, begin_table)
                            found_begin_table = true;
                        end
                    else
                        % Collect lines after finding the marker
                        % For example, you can save the next three lines as a single string
                        if or(or(strcmp(type,'H.2'), strcmp(type,'H.1')), strcmp(type,'H.0')) 
                            if line_counter < 3
                                if contains(line, 'flag')
                                    line = '';
                                end 
                                string = [string, line];
                                
                                if ~isempty(line)
                                    line_counter = line_counter + 1;
                                end
                            else
                                break; % Exit the loop after saving the next three lines
                            end
                        else
                            if line_counter<3*9
                                if ~isempty(line)
                                    l = strip(line);
                                    if ~isempty(str2num(l(1)))
                                        line_counter = line_counter + 1;
                                        string = [string,line];
                                    end
                                end
                            end
                        end
                    end
                end   
            end
        end
    end
    % Close the file
    fclose(fileID);

    string = strrep(string, 'D', 'E');
    numeric_values = regexp(string, '-?\d+\.\d+E[+-]\d+', 'match');

    if isempty(numeric_values)
        numeric_values = regexp(string, '-?\d+\.\d+e[+-]\d+', 'match');
    end
    
    if ~or(or(strcmp(type, 'H.1'), strcmp(type, 'H.2')), strcmp(type,'H.0'))
        num_arr = str2double(numeric_values);
        num_2 = zeros(27,3);
        
        for j=1:3
            for i=1:27
                num_2(i,j) = num_arr(j+3*(i-1));
            end
        end
        
        number_array = zeros(9);
        for i = 1:9
            for j=1:3
                number_array(i,j) = num_2(i,j);
                number_array(i,j+3) = num_2(i+9, j);
                number_array(i,j+6) = num_2(i+18, j);
            end
        end


    else 
        % Convert the extracted strings to double and store them in an array
        number_array = str2double(numeric_values);
    end


    
    
    if ~found_reaction
        number_array = zeros(1,9);
    end
end






