#!/usr/bin/env python3
"""
DEPRECATED: This file has been merged into app.py with hyperscan integration.

Use app.py instead which now includes:
- Modern Streamlit UI and styling
- Hyperscan-based motif detection via hyperscan_integration.py
- All functionality from this file

To run the application:
    streamlit run app.py

This file is kept for reference but should not be used.
"""

import streamlit as st

st.warning("⚠️ This Streamlit app has been deprecated!")
st.info("Please use `streamlit run app.py` instead.")
st.markdown("""
### Why was this deprecated?

This file (app_streamlit.py) has been successfully merged into app.py with the following improvements:

1. **Modern UI**: Uses the current style of Streamlit with professional CSS
2. **Hyperscan Integration**: Maintains the high-performance hyperscan approach
3. **Enhanced Features**: Comprehensive visualization and analysis tools
4. **Better Organization**: Cleaner code structure and better maintainability

### How to use the new app:

```bash
streamlit run app.py
```

The new app includes all the functionality from this file plus many enhancements.
""")