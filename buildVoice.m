function [noteNums,durs] = buildVoice(voiceIn,PulsesPerSong)
% Expand voice description.
% The input data only tells us the pulse positions where notes
% start and their durations. We use this to build a list of
% pulses for the whole piece. If no note starts at pulse
% k, then we set the note and duration for that pulse to zero.
durs = zeros(1,PulsesPerSong); % Pre-allocate array
noteNums = durs; % Pre-allocate array
for i=1:numel(voiceIn.startPulses)
    k = voiceIn.startPulses(i);
    noteNums(k) = voiceIn.noteNumbers(i);
    durs(k) = voiceIn.durations(i);
end