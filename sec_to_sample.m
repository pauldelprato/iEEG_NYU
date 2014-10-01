function sample =  sec_to_sample(time_s)
    Fs = 512;
    sample = fix(time_s*Fs)+1;
end