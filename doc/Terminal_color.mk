#Linux terminal color test
all_color: regular bold underline background
.PHONYP:all_color
regular:
	@echo -e "\e[0;30m Gray -Regular"
	@echo -e "\e[0;31m Red"
	@echo -e "\e[0;32m Green"
	@echo -e "\e[0;33m Yellow"
	@echo -e "\e[0;34m Blue"
	@echo -e "\e[0;35m Purple"
	@echo -e "\e[0;36m Cyan"
	@echo -e "\e[0;37m White"
	@echo -e "\e[0m Text Reset"
.PHONYP:regular
bold:
	@echo -e "\e[1;30m Gray -Bold"
	@echo -e "\e[1;31m Red"
	@echo -e "\e[1;32m Green"
	@echo -e "\e[1;33m Yellow"
	@echo -e "\e[1;34m Blue"
	@echo -e "\e[1;35m Purple"
	@echo -e "\e[1;36m Cyan"
	@echo -e "\e[1;37m White"
	@echo -e "\e[0m Text Reset"
.PHONYP:bold
underline:
	@echo -e "\e[4;30m Gray -Underline"
	@echo -e "\e[4;31m Red"
	@echo -e "\e[4;32m Green"
	@echo -e "\e[4;33m Yellow"
	@echo -e "\e[4;34m Blue"
	@echo -e "\e[4;35m Purple"
	@echo -e "\e[4;36m Cyan"
	@echo -e "\e[4;37m White"
	@echo -e "\e[0m Text Reset"
.PHONYP:underline
background:
	@echo -e "\e[40m Gray -Background"
	@echo -e "\e[41m Red"
	@echo -e "\e[42m Green"
	@echo -e "\e[43m Yellow"
	@echo -e "\e[44m Blue"
	@echo -e "\e[45m Purple"
	@echo -e "\e[46m Cyan"
	@echo -e "\e[47m White"
	@echo -e "\e[0m Text Reset"
.PHONYP:background
