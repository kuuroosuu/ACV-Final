function buttonH=wavePitchPlayButton(wObj, pitchObj1, pitchObj2)
% wavePitchPlayButton: Place two buttons for playing wave and pitch
%
%	Usage:
%		wavePitchPlayButton(wObj, pitchObj1, pitchObj2)
%
%	Description:
%		This function is used in waveAssess.m.

%	Roger Jang, 20051011, 20101007

if nargin<2, pitchObj1=[]; end
if nargin<3, pitchObj2=[]; end

% Normalized the amplitude if necessary
if ~wObj.amplitudeNormalized, wObj.signal=double(wObj.signal)/(2^wObj.nbits/2); wObj.amplitudeNormalized=1; end
% Assign to U for userData
U.wObj=wObj; U.pitch1=pitchObj1; U.pitch2=pitchObj2;

if ~isempty(wObj)
	buttonLabel='Play wave';
	if isfield(wObj, 'buttonLabel'), buttonLabel=wObj.buttonLabel; end
	buttonH(1)=uicontrol('string', 'Play Wave', 'position', [20, 10, 100, 20], 'userData', U);
	set(buttonH(1), 'Callback', 'U=get(gco, ''userData''); sound(U.wObj.signal, U.wObj.fs)');
%	set(buttonH(1), 'Callback', 'U=get(gco, ''userData''); audioPlayWithBar(U.wObj)');
end

if ~isempty(pitchObj1)
	buttonLabel='Play pitch';
	if isfield(pitchObj1, 'buttonLabel'), buttonLabel=pitchObj1.buttonLabel; end
	buttonH(2)=uicontrol('string', buttonLabel, 'position', [140, 10, 100, 20], 'userData', U);
%	set(buttonH(2), 'Callback', 'U=get(gco, ''userData''); xlim=get(gca, ''xlim''); pvPlay(U.pitch1.signal(find(xlim(1)*U.pitch1.frameRate<=1:length(U.pitch1.signal) & 1:length(U.pitch1.signal)<=xlim(2)*U.pitch1.frameRate)), U.pitch1.frameRate);');
	set(buttonH(2), 'Callback', 'U=get(gco, ''userData''); pvPlay(U.pitch1.signal, U.pitch1.frameRate);');
end

if ~isempty(wObj) & ~isempty(pitchObj1)
	buttonLabel='Play both';
	buttonH(3)=uicontrol('string', buttonLabel, 'position', [260, 10, 100, 20], 'userData', U);
	set(buttonH(3), 'Callback', 'U=get(gco, ''userData''); wavePitchPlay(U.wObj, U.pitch1);');
end

if ~isempty(pitchObj2)
	buttonLabel='Play Pitch2';
	if isfield(pitchObj2, 'buttonLabel'), buttonLabel=pitchObj2.buttonLabel; end
	buttonH(4)=uicontrol('string', buttonLabel, 'position', [380, 10, 100, 20], 'userData', U);
%	set(buttonH(4), 'Callback', 'U=get(gco, ''userData''); xlim=get(gca, ''xlim''); pvPlay(U.pitch2.signal(find(xlim(1)*U.pitch2.frameRate<=1:length(U.pitch2.signal) & 1:length(U.pitch2.signal)<=xlim(2)*U.pitch2.frameRate)), U.pitch2.frameRate);');
	set(buttonH(4), 'Callback', 'U=get(gco, ''userData''); pvPlay(U.pitch2.signal, U.pitch2.frameRate);');
end
