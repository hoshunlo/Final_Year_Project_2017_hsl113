function f = ClearMissingPupil(Pupil,N)
    temp_pupil = Pupil;
    for i = 1:N;
        if(Pupil(i) == -1)
            temp_pupil(i) = NaN;
        end
    end
    f = temp_pupil;
end
