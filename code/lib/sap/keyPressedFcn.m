function key_pressed_fcn(fig_obj, eventData)
%fprintf('CurrentKey = %s\n', get(fig_obj, 'CurrentKey'));
%fprintf('CurrentCharacter = %s\n', get(fig_obj, 'CurrentCharacter'));
%fprintf('CurrentModifier = %s\n\n', cell2str(get(fig_obj, 'CurrentModifier')));
%disp(eventData);
fprintf('Time=%s\n', get(findobj(gcf, 'tag', 'wavePlayButton'), 'string'));