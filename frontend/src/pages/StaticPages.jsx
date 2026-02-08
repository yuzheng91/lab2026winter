import React from 'react';
import { Container, Typography, Box } from '@mui/material';

export const Home = () => (
  <Container sx={{ mt: 5, textAlign: 'center' }}>
    <Typography variant="h3" gutterBottom>Welcome to WGLA</Typography>
    <Typography variant="h6" color="text.secondary">
      Please select "Tool" to start your gene analysis.
    </Typography>
  </Container>
);

export const Help = () => (
  <Container sx={{ mt: 5 }}>
    <Typography variant="h4">User Guide</Typography>
    <Typography sx={{ mt: 2 }}>Here you can place documentation...</Typography>
  </Container>
);

export const Contact = () => (
  <Container sx={{ mt: 5 }}>
    <Typography variant="h4">Contact Us</Typography>
    <Typography sx={{ mt: 2 }}>Email: support@ymla.org</Typography>
  </Container>
);