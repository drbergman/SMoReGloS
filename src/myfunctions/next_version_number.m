function num = next_version_number(formatSpec)

num = 1;
while exist(sprintf(formatSpec,num),'file')
    num = num+1;
end