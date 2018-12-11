durations = [1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.0, 6.0, 7.5, 9.0, 10.5, 12.0, 12.5, 15.0]
min_periods = [0.5, 0.52, 0.6517, 0.7824, 0.912, 1.178, 1.3056, 1.567, 1.952, 2.343, 2.75, 3.14, 3.257, 3.91]
#max_periods = [50.21, 118.147, 231.549, 397.82, 635.29
max_periods = [50.045, 118.626, 231.69, 400.359, 635.76, 725, 725, 725, 725, 725, 725, 725, 725, 725]
num_dur = 14

function get_legal_durations(period::Float64,duration::Float64)
  min_duration = 0.0
  max_duration = 0.0
#  for i in 1:num_durations
  i = 1
  while min_duration == 0.0 || max_duration == 0.0
#println(i," ",min_periods[num_dur+1-i])
    if period < max_periods[i] && min_duration == 0.0
      min_duration = durations[i]
    end 
    if period > min_periods[num_dur+1-i] && max_duration == 0.0
      max_duration = durations[num_dur+1-i]
    end
    if i > 14
      println("No durations match this period")
      return(0.0,0.0)
    end
    i+=1
  end
#  println(min_duration," ",max_duration)
  if duration<=max_duration && duration>=min_duration
    return true
  else
    return false
  end
#  return (min_duration,max_duration)
end
