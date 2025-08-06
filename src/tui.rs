use crossterm::{
    event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyCode, KeyEventKind},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use csv::ReaderBuilder;
use ratatui::{
    backend::{Backend, CrosstermBackend},
    layout::{Constraint, Direction, Layout},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Clear, List, ListItem, ListState, Paragraph},
    Frame, Terminal,
};
use std::{error::Error, fs, io, process::Command};
use toml_edit::{value, DocumentMut};

#[derive(Debug, Clone)]
pub struct PeptideEntry {
    pub name: String,         // PDB ID
    pub peptide_type: String, // Type (AH_MP, etc.)
    pub structured_region: String,
    pub sequence: String, // FASTA sequence
    pub total_length: u32,
    pub structured_length: u32,
    pub description: String,
    pub af2_min_rmsd: f64,
    pub af2_max_rmsd: f64,
    pub af2_avg_rmsd: f64,
    pub af2_std_rmsd: f64,
}

pub struct PeptideSelector {
    peptides: Vec<PeptideEntry>,
    filtered_indices: Vec<usize>,
    list_state: ListState,
    search_query: String,
    search_mode: bool,
    show_help: bool,
}

impl PeptideSelector {
    pub fn new() -> Result<Self, Box<dyn Error>> {
        let peptides = Self::load_peptides_from_csv("data.csv")?;
        let filtered_indices: Vec<usize> = (0..peptides.len()).collect();

        let mut list_state = ListState::default();
        if !filtered_indices.is_empty() {
            list_state.select(Some(0));
        }

        Ok(Self {
            peptides,
            filtered_indices,
            list_state,
            search_query: String::new(),
            search_mode: false,
            show_help: false,
        })
    }

    fn load_peptides_from_csv(file_path: &str) -> Result<Vec<PeptideEntry>, Box<dyn Error>> {
        let content = fs::read_to_string(file_path)?;
        let mut reader = ReaderBuilder::new()
            .has_headers(true)
            .from_reader(content.as_bytes());

        let mut peptides = Vec::new();

        for result in reader.records() {
            let record = result?;

            // Skip empty rows or rows with empty names
            if record.len() < 11 || record.get(1).unwrap_or("").trim().is_empty() {
                continue;
            }

            let name = record.get(1).unwrap_or("").trim().to_string();
            let peptide_type = record.get(2).unwrap_or("").trim().to_string();
            let structured_region = record.get(3).unwrap_or("").trim().to_string();
            let sequence = record.get(4).unwrap_or("").trim().to_string();
            let total_length = record.get(5).unwrap_or("0").trim().parse().unwrap_or(0);
            let structured_length = record.get(6).unwrap_or("0").trim().parse().unwrap_or(0);
            let description = record.get(7).unwrap_or("").trim().to_string();
            let af2_min_rmsd = record.get(8).unwrap_or("0.0").trim().parse().unwrap_or(0.0);
            let af2_max_rmsd = record.get(9).unwrap_or("0.0").trim().parse().unwrap_or(0.0);
            let af2_avg_rmsd = record
                .get(10)
                .unwrap_or("0.0")
                .trim()
                .parse()
                .unwrap_or(0.0);
            let af2_std_rmsd = record
                .get(11)
                .unwrap_or("0.0")
                .trim()
                .parse()
                .unwrap_or(0.0);

            // Skip entries with empty names or sequences
            if name.is_empty() || sequence.is_empty() {
                continue;
            }

            peptides.push(PeptideEntry {
                name,
                peptide_type,
                structured_region,
                sequence,
                total_length,
                structured_length,
                description,
                af2_min_rmsd,
                af2_max_rmsd,
                af2_avg_rmsd,
                af2_std_rmsd,
            });
        }

        Ok(peptides)
    }

    fn filter_peptides(&mut self) {
        if self.search_query.is_empty() {
            self.filtered_indices = (0..self.peptides.len()).collect();
        } else {
            let query = self.search_query.to_lowercase();
            self.filtered_indices = self
                .peptides
                .iter()
                .enumerate()
                .filter(|(_, peptide)| {
                    peptide.name.to_lowercase().contains(&query)
                        || peptide.sequence.to_lowercase().contains(&query)
                        || peptide.peptide_type.to_lowercase().contains(&query)
                        || peptide.description.to_lowercase().contains(&query)
                })
                .map(|(i, _)| i)
                .collect();
        }

        // Reset selection to first item
        if !self.filtered_indices.is_empty() {
            self.list_state.select(Some(0));
        } else {
            self.list_state.select(None);
        }
    }

    fn next(&mut self) {
        if self.filtered_indices.is_empty() {
            return;
        }

        let i = match self.list_state.selected() {
            Some(i) => {
                if i >= self.filtered_indices.len() - 1 {
                    0
                } else {
                    i + 1
                }
            }
            None => 0,
        };
        self.list_state.select(Some(i));
    }

    fn previous(&mut self) {
        if self.filtered_indices.is_empty() {
            return;
        }

        let i = match self.list_state.selected() {
            Some(i) => {
                if i == 0 {
                    self.filtered_indices.len() - 1
                } else {
                    i - 1
                }
            }
            None => 0,
        };
        self.list_state.select(Some(i));
    }

    fn get_selected_peptide(&self) -> Option<&PeptideEntry> {
        if let Some(selected_idx) = self.list_state.selected() {
            if let Some(&peptide_idx) = self.filtered_indices.get(selected_idx) {
                return self.peptides.get(peptide_idx);
            }
        }
        None
    }

    fn update_toml_and_run(&self) -> Result<(), Box<dyn Error>> {
        if let Some(peptide) = self.get_selected_peptide() {
            // Update TOML file
            let toml_content = fs::read_to_string("dncs.toml")?;
            let mut doc = toml_content.parse::<DocumentMut>()?;

            // Update only moleculename and sequence
            doc["simulation"]["moleculename"] = value(&peptide.name);
            doc["simulation"]["sequence"] = value(&peptide.sequence);

            // Write back to file
            fs::write("dncs.toml", doc.to_string())?;

            // Run the command (assuming it's a Python script)
            let output = Command::new("python").arg("python/src/main.py").output()?;

            if !output.status.success() {
                eprintln!(
                    "Command failed with error: {}",
                    String::from_utf8_lossy(&output.stderr)
                );
            } else {
                println!("Command executed successfully!");
                println!("Output: {}", String::from_utf8_lossy(&output.stdout));
            }
        }
        Ok(())
    }
}

impl PeptideSelector {
    pub fn run(&mut self) -> Result<(), Box<dyn Error>> {
        // Setup terminal
        enable_raw_mode()?;
        let mut stdout = io::stdout();
        execute!(stdout, EnterAlternateScreen, EnableMouseCapture)?;
        let backend = CrosstermBackend::new(stdout);
        let mut terminal = Terminal::new(backend)?;

        let result = self.run_app(&mut terminal);

        // Restore terminal
        disable_raw_mode()?;
        execute!(
            terminal.backend_mut(),
            LeaveAlternateScreen,
            DisableMouseCapture
        )?;
        terminal.show_cursor()?;

        result
    }

    fn run_app<B: Backend>(&mut self, terminal: &mut Terminal<B>) -> Result<(), Box<dyn Error>> {
        loop {
            terminal.draw(|f| self.ui(f))?;

            if let Event::Key(key) = event::read()? {
                if key.kind == KeyEventKind::Press {
                    match key.code {
                        KeyCode::Char('q') | KeyCode::Esc if !self.search_mode => {
                            return Ok(());
                        }
                        KeyCode::Char('h') | KeyCode::F(1) if !self.search_mode => {
                            self.show_help = !self.show_help;
                        }
                        KeyCode::Char('/') if !self.search_mode => {
                            self.search_mode = true;
                        }
                        KeyCode::Esc if self.search_mode => {
                            self.search_mode = false;
                            self.search_query.clear();
                            self.filter_peptides();
                        }
                        KeyCode::Enter if self.search_mode => {
                            self.search_mode = false;
                            self.filter_peptides();
                        }
                        KeyCode::Enter if !self.search_mode => {
                            // Update TOML and run command
                            if let Err(e) = self.update_toml_and_run() {
                                eprintln!("Error: {}", e);
                            }
                            return Ok(());
                        }
                        KeyCode::Down | KeyCode::Char('j') if !self.search_mode => {
                            self.next();
                        }
                        KeyCode::Up | KeyCode::Char('k') if !self.search_mode => {
                            self.previous();
                        }
                        KeyCode::Char(c) if self.search_mode => {
                            self.search_query.push(c);
                        }
                        KeyCode::Backspace if self.search_mode => {
                            self.search_query.pop();
                        }
                        _ => {}
                    }
                }
            }
        }
    }

    fn ui(&mut self, f: &mut Frame) {
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3), // Search bar
                Constraint::Min(0),    // Main list
                Constraint::Length(3), // Status bar
            ])
            .split(f.area());

        // Search bar
        let search_style = if self.search_mode {
            Style::default()
                .fg(Color::Yellow)
                .add_modifier(Modifier::BOLD)
        } else {
            Style::default().fg(Color::Gray)
        };

        let search_text = if self.search_mode {
            format!("Search: {}_", self.search_query)
        } else {
            format!("Search: {} (Press '/' to search)", self.search_query)
        };

        let search_paragraph = Paragraph::new(search_text).style(search_style).block(
            Block::default()
                .borders(Borders::ALL)
                .title("DNCS Peptide Selector"),
        );
        f.render_widget(search_paragraph, chunks[0]);

        // Main list
        let items: Vec<ListItem> = self
            .filtered_indices
            .iter()
            .map(|&i| {
                let peptide = &self.peptides[i];
                let content = format!(
                    "{:<8} | {:<30} | {:<8} | Len: {:<3} | {}",
                    peptide.name,
                    if peptide.sequence.len() > 30 {
                        format!("{}...", &peptide.sequence[..27])
                    } else {
                        peptide.sequence.clone()
                    },
                    peptide.peptide_type,
                    peptide.total_length,
                    if peptide.description.len() > 40 {
                        format!("{}...", &peptide.description[..37])
                    } else {
                        peptide.description.clone()
                    }
                );
                ListItem::new(content)
            })
            .collect();

        let list = List::new(items)
            .block(Block::default().borders(Borders::ALL).title(format!(
                "Peptides ({}/{})",
                self.filtered_indices.len(),
                self.peptides.len()
            )))
            .highlight_style(
                Style::default()
                    .bg(Color::Blue)
                    .fg(Color::White)
                    .add_modifier(Modifier::BOLD),
            )
            .highlight_symbol("> ");

        f.render_stateful_widget(&list, chunks[1], &mut self.list_state);

        // Status bar
        let status_text = if self.search_mode {
            "ESC: Exit search | ENTER: Apply filter"
        } else {
            "↑/↓: Navigate | ENTER: Select & Run | /: Search | h: Help | q: Quit"
        };

        let status = Paragraph::new(status_text)
            .style(Style::default().fg(Color::Cyan))
            .block(Block::default().borders(Borders::ALL));
        f.render_widget(status, chunks[2]);

        // Help popup
        if self.show_help {
            let help_text = vec![
                Line::from("DNCS Peptide Selector - Help"),
                Line::from(""),
                Line::from("Navigation:"),
                Line::from("  ↑/↓ or j/k  - Move up/down"),
                Line::from("  ENTER       - Select peptide and run"),
                Line::from(""),
                Line::from("Search:"),
                Line::from("  /           - Start search"),
                Line::from("  ESC         - Cancel search"),
                Line::from("  ENTER       - Apply search filter"),
                Line::from(""),
                Line::from("Other:"),
                Line::from("  h or F1     - Toggle this help"),
                Line::from("  q or ESC    - Quit application"),
                Line::from(""),
                Line::from("Press any key to close help..."),
            ];

            let help_paragraph = Paragraph::new(help_text)
                .block(Block::default().borders(Borders::ALL).title("Help"))
                .style(Style::default().bg(Color::Black).fg(Color::White));

            let area = centered_rect(60, 70, f.area());
            f.render_widget(Clear, area);
            f.render_widget(help_paragraph, area);
        }

        // Show selected peptide details in a side panel if there's space
        if f.area().width > 120 {
            let main_chunks = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([Constraint::Min(80), Constraint::Length(40)])
                .split(chunks[1]);

            // Re-render the list in the left panel
            f.render_stateful_widget(&list, main_chunks[0], &mut self.list_state);

            // Show details in the right panel
            if let Some(peptide) = self.get_selected_peptide() {
                let details = vec![
                    Line::from(vec![
                        Span::styled("Name: ", Style::default().add_modifier(Modifier::BOLD)),
                        Span::raw(&peptide.name),
                    ]),
                    Line::from(vec![
                        Span::styled("Type: ", Style::default().add_modifier(Modifier::BOLD)),
                        Span::raw(&peptide.peptide_type),
                    ]),
                    Line::from(vec![
                        Span::styled("Length: ", Style::default().add_modifier(Modifier::BOLD)),
                        Span::raw(peptide.total_length.to_string()),
                    ]),
                    Line::from(""),
                    Line::from(vec![Span::styled(
                        "Sequence:",
                        Style::default().add_modifier(Modifier::BOLD),
                    )]),
                    Line::from(Span::raw(&peptide.sequence)),
                    Line::from(""),
                    Line::from(vec![Span::styled(
                        "Description:",
                        Style::default().add_modifier(Modifier::BOLD),
                    )]),
                    Line::from(Span::raw(&peptide.description)),
                    Line::from(""),
                    Line::from(vec![Span::styled(
                        "RMSD Stats:",
                        Style::default().add_modifier(Modifier::BOLD),
                    )]),
                    Line::from(format!("  Min: {:.4}", peptide.af2_min_rmsd)),
                    Line::from(format!("  Max: {:.4}", peptide.af2_max_rmsd)),
                    Line::from(format!("  Avg: {:.4}", peptide.af2_avg_rmsd)),
                    Line::from(format!("  Std: {:.4}", peptide.af2_std_rmsd)),
                ];

                let details_paragraph = Paragraph::new(details)
                    .block(Block::default().borders(Borders::ALL).title("Details"))
                    .wrap(ratatui::widgets::Wrap { trim: true });
                f.render_widget(details_paragraph, main_chunks[1]);
            }
        }
    }
}

// Helper function to create a centered rectangle
fn centered_rect(
    percent_x: u16,
    percent_y: u16,
    r: ratatui::layout::Rect,
) -> ratatui::layout::Rect {
    let popup_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Percentage((100 - percent_y) / 2),
            Constraint::Percentage(percent_y),
            Constraint::Percentage((100 - percent_y) / 2),
        ])
        .split(r);

    Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage((100 - percent_x) / 2),
            Constraint::Percentage(percent_x),
            Constraint::Percentage((100 - percent_x) / 2),
        ])
        .split(popup_layout[1])[1]
}

pub fn run_peptide_selector() -> Result<(), Box<dyn Error>> {
    let mut selector = PeptideSelector::new()?;
    selector.run()
}
