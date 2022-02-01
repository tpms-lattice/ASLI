#ifndef GETLINE_MULTIOS_H
#define GETLINE_MULTIOS_H

inline std::istream& getline_multiOS(std::istream& is, std::string& s) {
	/* getline_multiOS is a getline substitute that handles line endings "\n", 
	*  "\r" and "\r\n".
	*/

	s.clear();

	std::istream::sentry se(is, true);
	std::streambuf* isbuf = is.rdbuf();

	while(true) {
		int ch = isbuf->sbumpc();
		switch (ch) {
			case '\n': // "\n" line endings
				return is;

			case '\r': // "\r" and "\r\n" line endings
				if(isbuf->sgetc() == '\n')
					isbuf->sbumpc();
				return is;

			case std::streambuf::traits_type::eof(): // End line without line ending
				if(s.empty())
					is.setstate(std::ios::eofbit);
				return is;

			default:
				s = s + (char)ch;
		}
	}
}

#endif