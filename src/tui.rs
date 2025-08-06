use crossterm::{
    event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyCode, KeyEventKind},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use csv::ReaderBuilder;
use ratatui::{
    backend::{Backend, CrosstermBackend},
    layout::{Alignment, Constraint, Direction, Layout},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Clear, List, ListItem, ListState, Paragraph, Tabs},
    Frame, Terminal,
};
use std::{error::Error, fs, io, process::Command};
use toml_edit::{value, DocumentMut};

// Ayu Dark Theme Colors
struct AyuTheme;

impl AyuTheme {
    // Background colors
    const BG_SECONDARY: Color = Color::Rgb(25, 30, 35); // #191E23 - Secondary background
    const BG_ACCENT: Color = Color::Rgb(35, 40, 45); // #23282D - Accent background

    // Foreground colors
    const FG_PRIMARY: Color = Color::Rgb(230, 237, 243); // #E6EDF3 - Primary text
    const FG_SECONDARY: Color = Color::Rgb(171, 178, 191); // #ABB2BF - Secondary text
    const FG_MUTED: Color = Color::Rgb(92, 99, 112); // #5C6370 - Muted text

    // Accent colors
    const ORANGE: Color = Color::Rgb(255, 160, 102); // #FFA066 - Orange accent
    const YELLOW: Color = Color::Rgb(255, 213, 128); // #FFD580 - Yellow accent
    const GREEN: Color = Color::Rgb(186, 230, 126); // #BAE67E - Green accent
    const BLUE: Color = Color::Rgb(115, 184, 255); // #73B8FF - Blue accent
    const PURPLE: Color = Color::Rgb(209, 154, 255); // #D19AFF - Purple accent
    const CYAN: Color = Color::Rgb(95, 219, 255); // #5FDBFF - Cyan accent
    const RED: Color = Color::Rgb(242, 151, 142); // #F2978E - Red accent

    // UI specific colors
    const BORDER: Color = Color::Rgb(45, 50, 55); // #2D3237 - Border color
    const SELECTION: Color = Color::Rgb(55, 65, 75); // #37414B - Selection background
}

#[derive(Debug, Clone)]
pub struct PeptideEntry {
    pub name: String,         // PDB ID
    pub peptide_type: String, // Type (AH_MP, etc.)
    pub sequence: String,     // FASTA sequence
    pub total_length: u32,
    pub description: String,
}

#[derive(Debug, Clone)]
pub struct SimulationConfig {
    pub moleculename: String,
    pub folder: String,
    pub sequence: String,
    pub interface: String,
    pub n_samples: u32,
    pub md_simulation: u32,
    pub temp: f64,
    pub forcefield: Vec<String>,
    pub device: String,
    pub solvent: u32,
    pub steps: u32,
    pub gamma: f64,
    pub dt: f64,
    pub md_steps: u32,
    pub method: String,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            moleculename: "6RRO".to_string(),
            folder: "Result".to_string(),
            sequence: "GFIVKRFKILV".to_string(),
            interface: "openmm".to_string(),
            n_samples: 2048,
            md_simulation: 10,
            temp: 300.0,
            forcefield: vec![
                "amber14-all.xml".to_string(),
                "amber14/tip3p.xml".to_string(),
            ],
            device: "CUDA".to_string(),
            solvent: 100,
            steps: 5000,
            gamma: 1.0,
            dt: 0.002,
            md_steps: 1000,
            method: "fold".to_string(),
        }
    }
}

#[derive(Debug, PartialEq)]
enum AppMode {
    PeptideSelection,
    Configuration,
    ConfigurationReview,
    Help,
}

#[derive(Debug, PartialEq, Clone)]
enum SortMode {
    None,
    NameAsc,
    NameDesc,
    LengthAsc,
    LengthDesc,
    TypeAsc,
    TypeDesc,
}

#[derive(Debug, PartialEq, Clone)]
enum ConfigField {
    Moleculename,
    Folder,
    Sequence,
    Interface,
    NSamples,
    MdSimulation,
    Temp,
    Forcefield,
    Device,
    Solvent,
    Steps,
    Gamma,
    Dt,
    MdSteps,
    Method,
}

pub struct PeptideSelector {
    peptides: Vec<PeptideEntry>,
    filtered_indices: Vec<usize>,
    list_state: ListState,
    search_query: String,
    search_mode: bool,
    app_mode: AppMode,
    config: SimulationConfig,
    config_list_state: ListState,
    config_fields: Vec<ConfigField>,
    editing_field: Option<ConfigField>,
    edit_buffer: String,
    tab_index: usize,
    sort_mode: SortMode,
    enum_selection_mode: bool,
    enum_options: Vec<String>,
    enum_selected_index: usize,
    enum_list_state: ListState,
}

impl PeptideSelector {
    pub fn new() -> Result<Self, Box<dyn Error>> {
        let peptides = Self::load_peptides_from_csv("data.csv")?;
        let filtered_indices: Vec<usize> = (0..peptides.len()).collect();

        let mut list_state = ListState::default();
        if !filtered_indices.is_empty() {
            list_state.select(Some(2)); // Skip headers
        }

        let config = Self::load_config_from_toml("dncs.toml").unwrap_or_default();

        let config_fields = vec![
            ConfigField::Moleculename,
            ConfigField::Folder,
            ConfigField::Sequence,
            ConfigField::Interface,
            ConfigField::NSamples,
            ConfigField::MdSimulation,
            ConfigField::Temp,
            ConfigField::Forcefield,
            ConfigField::Device,
            ConfigField::Solvent,
            ConfigField::Steps,
            ConfigField::Gamma,
            ConfigField::Dt,
            ConfigField::MdSteps,
            ConfigField::Method,
        ];

        let mut config_list_state = ListState::default();
        config_list_state.select(Some(0));

        Ok(Self {
            peptides,
            filtered_indices,
            list_state,
            search_query: String::new(),
            search_mode: false,
            app_mode: AppMode::PeptideSelection,
            config,
            config_list_state,
            config_fields,
            editing_field: None,
            edit_buffer: String::new(),
            tab_index: 0,
            sort_mode: SortMode::None,
            enum_selection_mode: false,
            enum_options: Vec::new(),
            enum_selected_index: 0,
            enum_list_state: ListState::default(),
        })
    }

    fn load_config_from_toml(file_path: &str) -> Result<SimulationConfig, Box<dyn Error>> {
        let content = fs::read_to_string(file_path)?;
        let doc = content.parse::<DocumentMut>()?;

        let sim = &doc["simulation"];
        Ok(SimulationConfig {
            moleculename: sim["moleculename"].as_str().unwrap_or("6RRO").to_string(),
            folder: sim["folder"].as_str().unwrap_or("Result").to_string(),
            sequence: sim["sequence"]
                .as_str()
                .unwrap_or("GFIVKRFKILV")
                .to_string(),
            interface: sim["interface"].as_str().unwrap_or("openmm").to_string(),
            n_samples: sim["n_samples"].as_integer().unwrap_or(2048) as u32,
            md_simulation: sim["md_simulation"].as_integer().unwrap_or(10) as u32,
            temp: sim["temp"].as_float().unwrap_or(300.0),
            forcefield: sim["forcefield"]
                .as_array()
                .map(|arr| {
                    arr.iter()
                        .filter_map(|v| v.as_str().map(|s| s.to_string()))
                        .collect()
                })
                .unwrap_or_else(|| {
                    vec![
                        "amber14-all.xml".to_string(),
                        "amber14/tip3p.xml".to_string(),
                    ]
                }),
            device: sim["device"].as_str().unwrap_or("CUDA").to_string(),
            solvent: sim["solvent"].as_integer().unwrap_or(100) as u32,
            steps: sim["steps"].as_integer().unwrap_or(5000) as u32,
            gamma: sim["gamma"].as_float().unwrap_or(1.0),
            dt: sim["dt"].as_float().unwrap_or(0.002),
            md_steps: sim["md_steps"].as_integer().unwrap_or(1000) as u32,
            method: sim["method"].as_str().unwrap_or("fold").to_string(),
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
            let sequence = record.get(4).unwrap_or("").trim().to_string();
            let total_length = record
                .get(5)
                .unwrap_or("0")
                .trim()
                .parse::<f64>()
                .unwrap_or(0.0) as u32;
            let description = record.get(7).unwrap_or("").trim().to_string();

            // Skip entries with empty names or sequences
            if name.is_empty() || sequence.is_empty() {
                continue;
            }

            peptides.push(PeptideEntry {
                name,
                peptide_type,
                sequence,
                total_length,
                description,
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
                    // Always include length in text search
                    let length_str = peptide.total_length.to_string();
                    let matches_text = peptide.name.to_lowercase().contains(&query)
                        || peptide.sequence.to_lowercase().contains(&query)
                        || peptide.peptide_type.to_lowercase().contains(&query)
                        || peptide.description.to_lowercase().contains(&query)
                        || length_str.contains(&query);

                    // Also check for exact length match if query is a number
                    let matches_length = if let Ok(length_query) = query.parse::<u32>() {
                        peptide.total_length == length_query
                    } else {
                        false
                    };

                    matches_text || matches_length
                })
                .map(|(i, _)| i)
                .collect();
        }

        // Apply sorting
        self.sort_peptides();

        // Reset selection to first item (after headers)
        if !self.filtered_indices.is_empty() {
            self.list_state.select(Some(2)); // Skip headers
        } else {
            self.list_state.select(None);
        }
    }

    fn sort_peptides(&mut self) {
        match self.sort_mode {
            SortMode::None => {
                // Keep original order
            }
            SortMode::NameAsc => {
                self.filtered_indices
                    .sort_by(|&a, &b| self.peptides[a].name.cmp(&self.peptides[b].name));
            }
            SortMode::NameDesc => {
                self.filtered_indices
                    .sort_by(|&a, &b| self.peptides[b].name.cmp(&self.peptides[a].name));
            }
            SortMode::LengthAsc => {
                self.filtered_indices.sort_by(|&a, &b| {
                    self.peptides[a]
                        .total_length
                        .cmp(&self.peptides[b].total_length)
                });
            }
            SortMode::LengthDesc => {
                self.filtered_indices.sort_by(|&a, &b| {
                    self.peptides[b]
                        .total_length
                        .cmp(&self.peptides[a].total_length)
                });
            }
            SortMode::TypeAsc => {
                self.filtered_indices.sort_by(|&a, &b| {
                    self.peptides[a]
                        .peptide_type
                        .cmp(&self.peptides[b].peptide_type)
                });
            }
            SortMode::TypeDesc => {
                self.filtered_indices.sort_by(|&a, &b| {
                    self.peptides[b]
                        .peptide_type
                        .cmp(&self.peptides[a].peptide_type)
                });
            }
        }
    }

    fn cycle_sort_mode(&mut self) {
        self.sort_mode = match self.sort_mode {
            SortMode::None => SortMode::LengthAsc,
            SortMode::LengthAsc => SortMode::LengthDesc,
            SortMode::LengthDesc => SortMode::NameAsc,
            SortMode::NameAsc => SortMode::NameDesc,
            SortMode::NameDesc => SortMode::TypeAsc,
            SortMode::TypeAsc => SortMode::TypeDesc,
            SortMode::TypeDesc => SortMode::None,
        };
        self.filter_peptides(); // Re-apply filtering with new sort
    }

    fn get_sort_indicator(&self) -> &str {
        match self.sort_mode {
            SortMode::None => "",
            SortMode::NameAsc => " [Name ↑]",
            SortMode::NameDesc => " [Name ↓]",
            SortMode::LengthAsc => " [Length ↑]",
            SortMode::LengthDesc => " [Length ↓]",
            SortMode::TypeAsc => " [Type ↑]",
            SortMode::TypeDesc => " [Type ↓]",
        }
    }

    fn next(&mut self) {
        if self.filtered_indices.is_empty() {
            return;
        }

        let i = match self.list_state.selected() {
            Some(i) => {
                // Skip header rows (2) and wrap around
                let max_index = self.filtered_indices.len() + 1; // +2 for headers -1 for 0-based
                if i >= max_index {
                    2 // Start after headers
                } else {
                    i + 1
                }
            }
            None => 2, // Start after headers
        };
        self.list_state.select(Some(i));
    }

    fn previous(&mut self) {
        if self.filtered_indices.is_empty() {
            return;
        }

        let i = match self.list_state.selected() {
            Some(i) => {
                // Skip header rows (2) and wrap around
                if i <= 2 {
                    self.filtered_indices.len() + 1 // Last actual item
                } else {
                    i - 1
                }
            }
            None => 2, // Start after headers
        };
        self.list_state.select(Some(i));
    }

    fn get_selected_peptide(&self) -> Option<&PeptideEntry> {
        if let Some(selected_idx) = self.list_state.selected() {
            // Adjust for header rows (subtract 2)
            if selected_idx >= 2 {
                let adjusted_idx = selected_idx - 2;
                if let Some(&peptide_idx) = self.filtered_indices.get(adjusted_idx) {
                    return self.peptides.get(peptide_idx);
                }
            }
        }
        None
    }

    fn update_toml_and_run(&self) -> Result<(), Box<dyn Error>> {
        // Save current configuration to TOML
        self.save_config_to_toml()?;

        // Run the command (try python3 first, then python)
        let python_cmd = if Command::new("python3").arg("--version").output().is_ok() {
            "python3"
        } else {
            "python"
        };

        let output = Command::new(python_cmd)
            .arg("python/src/main.py")
            .output()?;

        if !output.status.success() {
            eprintln!(
                "Command failed with error: {}",
                String::from_utf8_lossy(&output.stderr)
            );
        } else {
            println!("Command executed successfully!");
            println!("Output: {}", String::from_utf8_lossy(&output.stdout));
        }
        Ok(())
    }

    fn run_just_command(&self) -> Result<(), Box<dyn Error>> {
        // Save current configuration to TOML
        self.save_config_to_toml()?;

        // Run the just run command
        let output = Command::new("just").arg("run").output()?;

        if !output.status.success() {
            eprintln!(
                "Command failed with error: {}",
                String::from_utf8_lossy(&output.stderr)
            );
        } else {
            println!("Simulation started successfully!");
            println!("Output: {}", String::from_utf8_lossy(&output.stdout));
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
                    // Handle mode-specific keys first (including enum selection)
                    let handled = match self.app_mode {
                        AppMode::PeptideSelection => self.handle_peptide_keys(key.code)?,
                        AppMode::Configuration => self.handle_config_keys(key.code)?,
                        AppMode::ConfigurationReview => self.handle_review_keys(key.code)?,
                        AppMode::Help => false, // Help mode falls through to global keys
                    };

                    // If mode-specific handler didn't handle the key, try global keys
                    if !handled {
                        match key.code {
                            KeyCode::Char('q') => {
                                return Ok(());
                            }
                            KeyCode::Esc => {
                                // Only exit on ESC if not in search mode and no editing/enum selection active
                                if !self.search_mode
                                    && self.editing_field.is_none()
                                    && !self.enum_selection_mode
                                {
                                    return Ok(());
                                }
                            }
                            KeyCode::Char('h') | KeyCode::F(1) => {
                                self.app_mode = if self.app_mode == AppMode::Help {
                                    AppMode::PeptideSelection
                                } else {
                                    AppMode::Help
                                };
                            }
                            KeyCode::Tab => {
                                // Only allow tab switching from main modes, not review mode
                                if self.app_mode != AppMode::ConfigurationReview {
                                    self.tab_index = (self.tab_index + 1) % 2;
                                    self.app_mode = if self.tab_index == 0 {
                                        AppMode::PeptideSelection
                                    } else {
                                        AppMode::Configuration
                                    };
                                }
                            }
                            KeyCode::F(5) => {
                                // Save configuration
                                if let Err(e) = self.save_config_to_toml() {
                                    eprintln!("Error saving config: {}", e);
                                }
                            }
                            KeyCode::F(9) => {
                                // Run simulation
                                if let Err(e) = self.update_toml_and_run() {
                                    eprintln!("Error: {}", e);
                                }
                                return Ok(());
                            }
                            _ => {}
                        }
                    }
                }
            }
        }
    }

    fn handle_peptide_keys(&mut self, key: KeyCode) -> Result<bool, Box<dyn Error>> {
        match key {
            KeyCode::Char('/') if !self.search_mode => {
                self.search_mode = true;
            }
            KeyCode::Char('s') if !self.search_mode => {
                self.cycle_sort_mode();
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
                // Update config with selected peptide and go to review page
                if let Some(peptide) = self.get_selected_peptide() {
                    let name = peptide.name.clone();
                    let sequence = peptide.sequence.clone();
                    self.config.moleculename = name;
                    self.config.sequence = sequence;
                    self.app_mode = AppMode::ConfigurationReview;
                }
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
            _ => return Ok(false), // Key not handled
        }
        Ok(true) // Key was handled
    }

    fn handle_config_keys(&mut self, key: KeyCode) -> Result<bool, Box<dyn Error>> {
        // Handle enum selection mode first - highest priority
        if self.enum_selection_mode && self.editing_field.is_some() {
            match key {
                KeyCode::Down | KeyCode::Char('j') => {
                    if !self.enum_options.is_empty() {
                        let next_index = if self.enum_selected_index < self.enum_options.len() - 1 {
                            self.enum_selected_index + 1
                        } else {
                            0 // Wrap to beginning
                        };
                        self.enum_selected_index = next_index;
                        self.enum_list_state.select(Some(next_index));
                        self.edit_buffer = self.enum_options[next_index].clone();
                    }
                    return Ok(true);
                }
                KeyCode::Up | KeyCode::Char('k') => {
                    if !self.enum_options.is_empty() {
                        let prev_index = if self.enum_selected_index > 0 {
                            self.enum_selected_index - 1
                        } else {
                            self.enum_options.len() - 1 // Wrap to end
                        };
                        self.enum_selected_index = prev_index;
                        self.enum_list_state.select(Some(prev_index));
                        self.edit_buffer = self.enum_options[prev_index].clone();
                    }
                    return Ok(true);
                }
                KeyCode::Enter => {
                    // Save the selected enum value
                    if !self.edit_buffer.is_empty() {
                        self.save_edited_field();
                    }
                    // Clean up enum selection state
                    self.editing_field = None;
                    self.edit_buffer.clear();
                    self.enum_selection_mode = false;
                    self.enum_list_state = ListState::default();
                    self.enum_options.clear();
                    self.enum_selected_index = 0;
                    return Ok(true);
                }
                KeyCode::Esc => {
                    // Cancel enum selection without saving
                    self.editing_field = None;
                    self.edit_buffer.clear();
                    self.enum_selection_mode = false;
                    self.enum_list_state = ListState::default();
                    self.enum_options.clear();
                    self.enum_selected_index = 0;
                    return Ok(true);
                }
                _ => return Ok(false), // Let other keys fall through
            }
        }

        // Handle regular text editing mode
        if self.editing_field.is_some() && !self.enum_selection_mode {
            match key {
                KeyCode::Enter => {
                    // Save the edited field
                    self.save_edited_field();
                    self.editing_field = None;
                    self.edit_buffer.clear();
                    return Ok(true);
                }
                KeyCode::Esc => {
                    // Cancel editing
                    self.editing_field = None;
                    self.edit_buffer.clear();
                    self.enum_selection_mode = false;
                    self.enum_list_state = ListState::default();
                    self.enum_options.clear();
                    return Ok(true);
                }
                KeyCode::Char(c) => {
                    self.edit_buffer.push(c);
                    return Ok(true);
                }
                KeyCode::Backspace => {
                    self.edit_buffer.pop();
                    return Ok(true);
                }
                _ => return Ok(false),
            }
        }

        // Regular navigation and interaction
        match key {
            KeyCode::Down | KeyCode::Char('j') => {
                self.next_config_field();
                Ok(true)
            }
            KeyCode::Up | KeyCode::Char('k') => {
                self.previous_config_field();
                Ok(true)
            }
            KeyCode::Enter => {
                // Start editing the selected field
                self.start_editing_field();
                Ok(true)
            }
            _ => Ok(false), // Key not handled
        }
    }

    fn handle_review_keys(&mut self, key: KeyCode) -> Result<bool, Box<dyn Error>> {
        // Handle enum selection mode first
        if self.enum_selection_mode && self.editing_field.is_some() {
            match key {
                KeyCode::Down | KeyCode::Char('j') => {
                    let next_index = if self.enum_selected_index < self.enum_options.len() - 1 {
                        self.enum_selected_index + 1
                    } else {
                        0 // Wrap to beginning
                    };
                    self.enum_selected_index = next_index;
                    self.enum_list_state.select(Some(next_index));
                    self.edit_buffer = self.enum_options[next_index].clone();
                }
                KeyCode::Up | KeyCode::Char('k') => {
                    let prev_index = if self.enum_selected_index > 0 {
                        self.enum_selected_index - 1
                    } else {
                        self.enum_options.len() - 1 // Wrap to end
                    };
                    self.enum_selected_index = prev_index;
                    self.enum_list_state.select(Some(prev_index));
                    self.edit_buffer = self.enum_options[prev_index].clone();
                }
                KeyCode::Enter => {
                    // Save the selected enum value
                    self.save_edited_field();
                    self.editing_field = None;
                    self.edit_buffer.clear();
                    self.enum_selection_mode = false;
                    self.enum_list_state = ListState::default();
                    self.enum_options.clear();
                }
                // Remove ESC exit behavior - user doesn't want it
                _ => {}
            }
            return Ok(true);
        }

        // Regular key handling
        match key {
            KeyCode::Down | KeyCode::Char('j') => {
                if self.editing_field.is_none() {
                    self.next_config_field();
                }
            }
            KeyCode::Up | KeyCode::Char('k') => {
                if self.editing_field.is_none() {
                    self.previous_config_field();
                }
            }
            KeyCode::Enter => {
                if self.editing_field.is_some() && !self.enum_selection_mode {
                    // Save the edited field
                    self.save_edited_field();
                    self.editing_field = None;
                    self.edit_buffer.clear();
                } else {
                    // Start editing the selected field
                    self.start_editing_field();
                }
            }
            KeyCode::Char('r') if self.editing_field.is_none() => {
                // Run simulation with current configuration
                if let Err(e) = self.run_just_command() {
                    eprintln!("Error running simulation: {}", e);
                }
                return Ok(true);
            }
            KeyCode::Esc => {
                if self.editing_field.is_some() {
                    // Cancel editing
                    self.editing_field = None;
                    self.edit_buffer.clear();
                    self.enum_selection_mode = false;
                    self.enum_list_state = ListState::default();
                    self.enum_options.clear();
                } else {
                    // Go back to peptide selection
                    self.app_mode = AppMode::PeptideSelection;
                }
            }
            KeyCode::Char(c) if self.editing_field.is_some() && !self.enum_selection_mode => {
                self.edit_buffer.push(c);
            }
            KeyCode::Backspace if self.editing_field.is_some() && !self.enum_selection_mode => {
                self.edit_buffer.pop();
            }
            _ => return Ok(false), // Key not handled
        }
        Ok(true) // Key was handled
    }

    fn next_config_field(&mut self) {
        let i = match self.config_list_state.selected() {
            Some(i) => {
                if i >= self.config_fields.len() - 1 {
                    0
                } else {
                    i + 1
                }
            }
            None => 0,
        };
        self.config_list_state.select(Some(i));
    }

    fn previous_config_field(&mut self) {
        let i = match self.config_list_state.selected() {
            Some(i) => {
                if i == 0 {
                    self.config_fields.len() - 1
                } else {
                    i - 1
                }
            }
            None => 0,
        };
        self.config_list_state.select(Some(i));
    }

    fn start_editing_field(&mut self) {
        if let Some(selected) = self.config_list_state.selected() {
            if let Some(field) = self.config_fields.get(selected) {
                self.editing_field = Some(field.clone());

                // Check if this field has enum options
                if let Some(options) = self.get_enum_options(field) {
                    self.enum_selection_mode = true;
                    self.enum_options = options;
                    // Find current value index
                    let current_value = self.get_field_value(field);
                    self.enum_selected_index = self
                        .enum_options
                        .iter()
                        .position(|opt| opt == &current_value)
                        .unwrap_or(0);
                    // Ensure the enum list state is properly initialized
                    let mut enum_state = ListState::default();
                    enum_state.select(Some(self.enum_selected_index));
                    self.enum_list_state = enum_state;
                    self.edit_buffer = self.enum_options[self.enum_selected_index].clone();
                } else {
                    self.enum_selection_mode = false;
                    self.enum_list_state = ListState::default();
                    self.edit_buffer = self.get_field_value(field);
                }
            }
        }
    }

    fn get_enum_options(&self, field: &ConfigField) -> Option<Vec<String>> {
        match field {
            ConfigField::Interface => Some(vec![
                "openmm".to_string(),
                "gromacs".to_string(),
                "amber".to_string(),
            ]),
            ConfigField::Device => Some(vec![
                "CPU".to_string(),
                "CUDA".to_string(),
                "HIP".to_string(),
                "OpenCL".to_string(),
            ]),
            ConfigField::Method => Some(vec![
                "fold".to_string(),
                "search".to_string(),
                "explore".to_string(),
            ]),
            ConfigField::Forcefield => Some(self.get_available_forcefields()),
            _ => None,
        }
    }

    fn get_available_forcefields(&self) -> Vec<String> {
        let mut forcefields = vec![
            "amber14-all.xml".to_string(),
            "amber14/tip3p.xml".to_string(),
            "amber14/tip4pew.xml".to_string(),
            "amber14/spce.xml".to_string(),
            "amber99sb.xml".to_string(),
            "amber99sbildn.xml".to_string(),
            "charmm36.xml".to_string(),
        ];

        // Try to scan library directory for more
        if let Ok(entries) = fs::read_dir("library/ForceFields") {
            for entry in entries.flatten() {
                if let Some(name) = entry.file_name().to_str() {
                    if name.ends_with(".xml") && !forcefields.contains(&name.to_string()) {
                        forcefields.push(name.to_string());
                    }
                }
            }
        }

        forcefields.sort();
        forcefields
    }

    fn check_molecule_exists(&self, molecule_name: &str) -> bool {
        // Check if the molecule directory exists in the results folder
        let result_path = format!("{}/{}", self.config.folder, molecule_name);
        std::path::Path::new(&result_path).exists()
    }

    fn get_molecule_status_indicator(&self, molecule_name: &str) -> (&str, Style) {
        if self.check_molecule_exists(molecule_name) {
            (
                "✓",
                Style::default()
                    .fg(AyuTheme::GREEN)
                    .add_modifier(Modifier::BOLD),
            )
        } else {
            ("○", Style::default().fg(AyuTheme::FG_MUTED))
        }
    }

    fn get_field_value(&self, field: &ConfigField) -> String {
        match field {
            ConfigField::Moleculename => self.config.moleculename.clone(),
            ConfigField::Folder => self.config.folder.clone(),
            ConfigField::Sequence => self.config.sequence.clone(),
            ConfigField::Interface => self.config.interface.clone(),
            ConfigField::NSamples => self.config.n_samples.to_string(),
            ConfigField::MdSimulation => self.config.md_simulation.to_string(),
            ConfigField::Temp => self.config.temp.to_string(),
            ConfigField::Forcefield => self.config.forcefield.join(", "),
            ConfigField::Device => self.config.device.clone(),
            ConfigField::Solvent => self.config.solvent.to_string(),
            ConfigField::Steps => self.config.steps.to_string(),
            ConfigField::Gamma => self.config.gamma.to_string(),
            ConfigField::Dt => self.config.dt.to_string(),
            ConfigField::MdSteps => self.config.md_steps.to_string(),
            ConfigField::Method => self.config.method.clone(),
        }
    }

    fn save_edited_field(&mut self) {
        if let Some(field) = &self.editing_field {
            match field {
                ConfigField::Moleculename => self.config.moleculename = self.edit_buffer.clone(),
                ConfigField::Folder => self.config.folder = self.edit_buffer.clone(),
                ConfigField::Sequence => self.config.sequence = self.edit_buffer.clone(),
                ConfigField::Interface => self.config.interface = self.edit_buffer.clone(),
                ConfigField::NSamples => {
                    if let Ok(val) = self.edit_buffer.parse() {
                        self.config.n_samples = val;
                    }
                }
                ConfigField::MdSimulation => {
                    if let Ok(val) = self.edit_buffer.parse() {
                        self.config.md_simulation = val;
                    }
                }
                ConfigField::Temp => {
                    if let Ok(val) = self.edit_buffer.parse() {
                        self.config.temp = val;
                    }
                }
                ConfigField::Forcefield => {
                    if self.enum_selection_mode {
                        // For enum mode, replace the entire forcefield list with the selected option
                        self.config.forcefield = vec![self.edit_buffer.clone()];
                    } else {
                        // For manual editing, allow comma-separated values
                        self.config.forcefield = self
                            .edit_buffer
                            .split(',')
                            .map(|s| s.trim().to_string())
                            .filter(|s| !s.is_empty())
                            .collect();
                    }
                }
                ConfigField::Device => self.config.device = self.edit_buffer.clone(),
                ConfigField::Solvent => {
                    if let Ok(val) = self.edit_buffer.parse() {
                        self.config.solvent = val;
                    }
                }
                ConfigField::Steps => {
                    if let Ok(val) = self.edit_buffer.parse() {
                        self.config.steps = val;
                    }
                }
                ConfigField::Gamma => {
                    if let Ok(val) = self.edit_buffer.parse() {
                        self.config.gamma = val;
                    }
                }
                ConfigField::Dt => {
                    if let Ok(val) = self.edit_buffer.parse() {
                        self.config.dt = val;
                    }
                }
                ConfigField::MdSteps => {
                    if let Ok(val) = self.edit_buffer.parse() {
                        self.config.md_steps = val;
                    }
                }
                ConfigField::Method => self.config.method = self.edit_buffer.clone(),
            }
        }
    }

    fn save_config_to_toml(&self) -> Result<(), Box<dyn Error>> {
        let toml_content = fs::read_to_string("dncs.toml")?;
        let mut doc = toml_content.parse::<DocumentMut>()?;

        doc["simulation"]["moleculename"] = value(&self.config.moleculename);
        doc["simulation"]["folder"] = value(&self.config.folder);
        doc["simulation"]["sequence"] = value(&self.config.sequence);
        doc["simulation"]["interface"] = value(&self.config.interface);
        doc["simulation"]["n_samples"] = value(self.config.n_samples as i64);
        doc["simulation"]["md_simulation"] = value(self.config.md_simulation as i64);
        doc["simulation"]["sample_top_n"] = value(self.config.md_simulation as i64);
        doc["simulation"]["temp"] = value(self.config.temp);
        doc["simulation"]["device"] = value(&self.config.device);
        doc["simulation"]["solvent"] = value(self.config.solvent as i64);
        doc["simulation"]["steps"] = value(self.config.steps as i64);
        doc["simulation"]["gamma"] = value(self.config.gamma);
        doc["simulation"]["dt"] = value(self.config.dt);
        doc["simulation"]["md_steps"] = value(self.config.md_steps as i64);
        doc["simulation"]["method"] = value(&self.config.method);

        // Handle forcefield array
        let mut ff_array = toml_edit::Array::new();
        for ff in &self.config.forcefield {
            ff_array.push(ff.as_str());
        }
        doc["simulation"]["forcefield"] = toml_edit::value(ff_array);

        fs::write("dncs.toml", doc.to_string())?;
        Ok(())
    }

    fn ui(&mut self, f: &mut Frame) {
        match self.app_mode {
            AppMode::Help => {
                self.render_help(f);
                return;
            }
            _ => {}
        }

        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3), // Tabs
                Constraint::Min(0),    // Main content
                Constraint::Length(3), // Status bar
            ])
            .split(f.area());

        // Render tabs
        let tab_titles = vec!["Peptide Selection", "Configuration"];
        let tabs = Tabs::new(tab_titles)
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .title("DNCS Control Panel")
                    .border_style(Style::default().fg(AyuTheme::BORDER))
                    .title_style(
                        Style::default()
                            .fg(AyuTheme::ORANGE)
                            .add_modifier(Modifier::BOLD),
                    ),
            )
            .style(Style::default().fg(AyuTheme::FG_SECONDARY))
            .highlight_style(
                Style::default()
                    .fg(AyuTheme::CYAN)
                    .bg(AyuTheme::SELECTION)
                    .add_modifier(Modifier::BOLD),
            )
            .select(self.tab_index);
        f.render_widget(tabs, chunks[0]);

        // Render main content based on current mode
        match self.app_mode {
            AppMode::PeptideSelection => self.render_peptide_selection(f, chunks[1]),
            AppMode::Configuration => self.render_configuration(f, chunks[1]),
            AppMode::ConfigurationReview => self.render_configuration_review(f, chunks[1]),
            AppMode::Help => {} // Already handled above
        }

        // Status bar
        let status_text = match self.app_mode {
            AppMode::PeptideSelection => {
                if self.search_mode {
                    "ESC: Exit search | ENTER: Apply filter"
                } else {
                    "↑/↓: Navigate | ENTER: Review & Run | /: Search | s: Sort | TAB: Config | F9: Run | h: Help | q: Quit"
                }
            }
            AppMode::Configuration => {
                if self.enum_selection_mode {
                    "▼ DROPDOWN ACTIVE: ↑/↓ Navigate | ENTER: Select | ESC: Cancel"
                } else if self.editing_field.is_some() {
                    "Type value | ENTER: Save | ESC: Cancel"
                } else {
                    "↑/↓: Navigate | ENTER: Edit | TAB: Peptides | F5: Save | F9: Run | h: Help | q: Quit"
                }
            }
            AppMode::ConfigurationReview => {
                if self.enum_selection_mode {
                    "↑/↓: Select option | ENTER: Confirm | ESC: Cancel"
                } else if self.editing_field.is_some() {
                    "Type value | ENTER: Save | ESC: Cancel editing"
                } else {
                    "↑/↓: Navigate | ENTER: Edit | r: RUN SIMULATION | ESC: Back | h: Help"
                }
            }
            AppMode::Help => "Any key to return",
        };

        let status = Paragraph::new(status_text)
            .style(Style::default().fg(AyuTheme::FG_PRIMARY))
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .border_style(Style::default().fg(AyuTheme::BORDER))
                    .title("Controls")
                    .title_style(Style::default().fg(AyuTheme::PURPLE)),
            );
        f.render_widget(status, chunks[2]);
    }

    fn render_peptide_selection(&mut self, f: &mut Frame, area: ratatui::layout::Rect) {
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3), // Search bar
                Constraint::Min(0),    // Main list
            ])
            .split(area);

        // Search bar
        let search_style = if self.search_mode {
            Style::default()
                .fg(AyuTheme::YELLOW)
                .bg(AyuTheme::BG_ACCENT)
                .add_modifier(Modifier::BOLD)
        } else {
            Style::default().fg(AyuTheme::FG_SECONDARY)
        };

        let search_text = if self.search_mode {
            format!("Search: {}_", self.search_query)
        } else {
            format!("Search: {} (Press '/' to search)", self.search_query)
        };

        let search_title = if self.search_mode {
            "Search Peptides"
        } else {
            "Search Peptides (✓ = Results Available, ○ = Not Processed)"
        };

        let search_paragraph = Paragraph::new(search_text).style(search_style).block(
            Block::default()
                .borders(Borders::ALL)
                .title(search_title)
                .border_style(Style::default().fg(if self.search_mode {
                    AyuTheme::YELLOW
                } else {
                    AyuTheme::BORDER
                }))
                .title_style(
                    Style::default()
                        .fg(AyuTheme::GREEN)
                        .add_modifier(Modifier::BOLD),
                ),
        );
        f.render_widget(search_paragraph, chunks[0]);

        // Main list - collect data first to avoid borrowing issues
        let peptide_data: Vec<_> = self
            .filtered_indices
            .iter()
            .map(|&i| {
                let peptide = &self.peptides[i];
                (
                    peptide.name.clone(),
                    peptide.sequence.clone(),
                    peptide.peptide_type.clone(),
                    peptide.total_length,
                    peptide.description.clone(),
                )
            })
            .collect();

        let items: Vec<ListItem> = peptide_data
            .iter()
            .map(
                |(name, sequence, peptide_type, total_length, description)| {
                    let (status_icon, status_style) = if self.check_molecule_exists(name) {
                        (
                            "✓",
                            Style::default()
                                .fg(AyuTheme::GREEN)
                                .add_modifier(Modifier::BOLD),
                        )
                    } else {
                        ("○", Style::default().fg(AyuTheme::FG_MUTED))
                    };

                    let formatted_name = format!("{:<6}", name);
                    let formatted_sequence = if sequence.len() > 30 {
                        format!("{:<30}", format!("{}...", &sequence[..27]))
                    } else {
                        format!("{:<30}", sequence)
                    };
                    let formatted_type = format!("{:<8}", peptide_type);
                    let formatted_length = format!("{:>3}", total_length);
                    let formatted_desc = if description.len() > 40 {
                        format!("{}...", &description[..37])
                    } else {
                        description.clone()
                    };

                    let content = vec![Line::from(vec![
                        Span::styled(status_icon, status_style),
                        Span::raw(" "),
                        Span::styled(
                            formatted_name,
                            Style::default()
                                .fg(AyuTheme::YELLOW)
                                .add_modifier(Modifier::BOLD),
                        ),
                        Span::raw(" │ "),
                        Span::styled(
                            formatted_sequence,
                            Style::default().fg(AyuTheme::FG_PRIMARY),
                        ),
                        Span::raw(" │ "),
                        Span::styled(formatted_type, Style::default().fg(AyuTheme::PURPLE)),
                        Span::raw(" │ Len:"),
                        Span::styled(formatted_length, Style::default().fg(AyuTheme::CYAN)),
                        Span::raw(" │ "),
                        Span::styled(formatted_desc, Style::default().fg(AyuTheme::FG_SECONDARY)),
                    ])];
                    ListItem::new(content)
                },
            )
            .collect();

        // Add header row for column alignment
        let mut items_with_header = Vec::new();

        // Create header
        let header_content = vec![Line::from(vec![
            Span::styled("S ", Style::default().fg(AyuTheme::FG_MUTED)),
            Span::styled(
                format!("{:<6}", "ID"),
                Style::default()
                    .fg(AyuTheme::FG_MUTED)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::styled(" │ ", Style::default().fg(AyuTheme::FG_MUTED)),
            Span::styled(
                format!("{:<30}", "SEQUENCE"),
                Style::default()
                    .fg(AyuTheme::FG_MUTED)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::styled(" │ ", Style::default().fg(AyuTheme::FG_MUTED)),
            Span::styled(
                format!("{:<8}", "TYPE"),
                Style::default()
                    .fg(AyuTheme::FG_MUTED)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::styled(" │ ", Style::default().fg(AyuTheme::FG_MUTED)),
            Span::styled(
                "LEN",
                Style::default()
                    .fg(AyuTheme::FG_MUTED)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::styled(" │ ", Style::default().fg(AyuTheme::FG_MUTED)),
            Span::styled(
                "DESCRIPTION",
                Style::default()
                    .fg(AyuTheme::FG_MUTED)
                    .add_modifier(Modifier::BOLD),
            ),
        ])];
        items_with_header
            .push(ListItem::new(header_content).style(Style::default().bg(AyuTheme::BG_ACCENT)));

        // Add separator line
        let separator_content = vec![Line::from("─".repeat(120))];
        items_with_header
            .push(ListItem::new(separator_content).style(Style::default().fg(AyuTheme::FG_MUTED)));

        // Add all peptide items
        items_with_header.extend(items);

        let list = List::new(items_with_header)
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .title(format!(
                        "Peptides ({}/{}){}",
                        self.filtered_indices.len(),
                        self.peptides.len(),
                        self.get_sort_indicator()
                    ))
                    .border_style(Style::default().fg(AyuTheme::BORDER))
                    .title_style(
                        Style::default()
                            .fg(AyuTheme::BLUE)
                            .add_modifier(Modifier::BOLD),
                    ),
            )
            .style(Style::default().fg(AyuTheme::FG_PRIMARY))
            .highlight_style(
                Style::default()
                    .bg(AyuTheme::SELECTION)
                    .fg(AyuTheme::GREEN)
                    .add_modifier(Modifier::BOLD),
            )
            .highlight_symbol("▶ ");
        f.render_stateful_widget(&list, chunks[1], &mut self.list_state);

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
                let (status_icon, status_style) = self.get_molecule_status_indicator(&peptide.name);
                let status_text = if self.check_molecule_exists(&peptide.name) {
                    "Results Available"
                } else {
                    "Not Processed"
                };

                let details = vec![
                    Line::from(vec![
                        Span::styled(
                            "Name: ",
                            Style::default()
                                .fg(AyuTheme::ORANGE)
                                .add_modifier(Modifier::BOLD),
                        ),
                        Span::styled(&peptide.name, Style::default().fg(AyuTheme::YELLOW)),
                    ]),
                    Line::from(vec![
                        Span::styled(
                            "Status: ",
                            Style::default()
                                .fg(AyuTheme::ORANGE)
                                .add_modifier(Modifier::BOLD),
                        ),
                        Span::styled(status_icon, status_style),
                        Span::raw(" "),
                        Span::styled(status_text, status_style),
                    ]),
                    Line::from(vec![
                        Span::styled(
                            "Type: ",
                            Style::default()
                                .fg(AyuTheme::ORANGE)
                                .add_modifier(Modifier::BOLD),
                        ),
                        Span::styled(&peptide.peptide_type, Style::default().fg(AyuTheme::PURPLE)),
                    ]),
                    Line::from(vec![
                        Span::styled(
                            "Length: ",
                            Style::default()
                                .fg(AyuTheme::ORANGE)
                                .add_modifier(Modifier::BOLD),
                        ),
                        Span::styled(
                            peptide.total_length.to_string(),
                            Style::default().fg(AyuTheme::CYAN),
                        ),
                    ]),
                    Line::from(""),
                    Line::from(vec![Span::styled(
                        "Sequence:",
                        Style::default()
                            .fg(AyuTheme::GREEN)
                            .add_modifier(Modifier::BOLD),
                    )]),
                    Line::from(Span::styled(
                        &peptide.sequence,
                        Style::default().fg(AyuTheme::FG_PRIMARY),
                    )),
                    Line::from(""),
                    Line::from(vec![Span::styled(
                        "Description:",
                        Style::default()
                            .fg(AyuTheme::GREEN)
                            .add_modifier(Modifier::BOLD),
                    )]),
                    Line::from(Span::styled(
                        &peptide.description,
                        Style::default().fg(AyuTheme::FG_SECONDARY),
                    )),
                ];

                let details_paragraph = Paragraph::new(details)
                    .block(
                        Block::default()
                            .borders(Borders::ALL)
                            .title("Details")
                            .border_style(Style::default().fg(AyuTheme::BORDER))
                            .title_style(
                                Style::default()
                                    .fg(AyuTheme::PURPLE)
                                    .add_modifier(Modifier::BOLD),
                            ),
                    )
                    .wrap(ratatui::widgets::Wrap { trim: true });
                f.render_widget(details_paragraph, main_chunks[1]);
            }
        }
    }

    fn render_configuration(&mut self, f: &mut Frame, area: ratatui::layout::Rect) {
        let chunks = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Percentage(50), Constraint::Percentage(50)])
            .split(area);

        // Use fixed width for consistent alignment - longest field is "MD Simulations (Top N)" = 22 chars

        // Configuration fields list
        let config_items: Vec<ListItem> = self
            .config_fields
            .iter()
            .map(|field| {
                let (name, value) = match field {
                    ConfigField::Moleculename => {
                        ("Molecule Name", self.config.moleculename.clone())
                    }
                    ConfigField::Folder => ("Output Folder", self.config.folder.clone()),
                    ConfigField::Sequence => (
                        "Sequence",
                        if self.config.sequence.len() > 20 {
                            format!("{}...", &self.config.sequence[..17])
                        } else {
                            self.config.sequence.clone()
                        },
                    ),
                    ConfigField::Interface => ("Interface ⚡", self.config.interface.clone()),
                    ConfigField::NSamples => {
                        ("Number of Samples", self.config.n_samples.to_string())
                    }
                    ConfigField::MdSimulation => (
                        "MD Simulations (Top N)",
                        self.config.md_simulation.to_string(),
                    ),
                    ConfigField::Temp => ("Temperature (K)", self.config.temp.to_string()),
                    ConfigField::Forcefield => (
                        "Force Fields ⚡",
                        if self.config.forcefield.len() <= 2 {
                            self.config.forcefield.join(", ")
                        } else {
                            format!(
                                "{} fields: {}, ...",
                                self.config.forcefield.len(),
                                &self.config.forcefield[0]
                            )
                        },
                    ),
                    ConfigField::Device => ("Device ⚡", self.config.device.clone()),
                    ConfigField::Solvent => ("Solvent Count", self.config.solvent.to_string()),
                    ConfigField::Steps => ("Equilibration Steps", self.config.steps.to_string()),
                    ConfigField::Gamma => ("Friction Coefficient", self.config.gamma.to_string()),
                    ConfigField::Dt => ("Time Step (ps)", self.config.dt.to_string()),
                    ConfigField::MdSteps => ("MD Steps", self.config.md_steps.to_string()),
                    ConfigField::Method => ("Sampling Method ⚡", self.config.method.clone()),
                };

                // Check if this field has enum options
                let has_enum_options = self.get_enum_options(field).is_some();

                let content = if Some(field) == self.editing_field.as_ref() {
                    if self.enum_selection_mode && has_enum_options {
                        // Show active dropdown selection with highlight
                        format!("{:<25} │ ▼ {} ◄ SELECTING", name, self.edit_buffer)
                    } else if has_enum_options {
                        // Show dropdown field ready for selection
                        format!("{:<25} │ ▼ {} ← PRESS ENTER", name, value)
                    } else {
                        // Regular text input field with cursor
                        format!("{:<25} │   {}_", name, self.edit_buffer)
                    }
                } else {
                    if has_enum_options {
                        // Show dropdown indicator for available enum fields
                        format!("{:<25} │ ▼ {}", name, value)
                    } else {
                        // Regular field display
                        format!("{:<25} │   {}", name, value)
                    }
                };

                ListItem::new(content)
            })
            .collect();

        let config_list = List::new(config_items)
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .title("Configuration")
                    .border_style(Style::default().fg(AyuTheme::BORDER))
                    .title_style(
                        Style::default()
                            .fg(AyuTheme::ORANGE)
                            .add_modifier(Modifier::BOLD),
                    ),
            )
            .style(Style::default().fg(AyuTheme::FG_PRIMARY))
            .highlight_style(
                Style::default()
                    .bg(AyuTheme::SELECTION)
                    .fg(AyuTheme::GREEN)
                    .add_modifier(Modifier::BOLD),
            )
            .highlight_symbol("▶ ");

        f.render_stateful_widget(&config_list, chunks[0], &mut self.config_list_state);

        // Render enum selection popup if in enum mode
        if self.enum_selection_mode && !self.enum_options.is_empty() {
            self.render_enum_popup(f, chunks[0]);
        }

        // Configuration details and help
        let details = vec![
            Line::from(vec![Span::styled(
                "Configuration Details",
                Style::default()
                    .fg(AyuTheme::CYAN)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Current Settings:",
                Style::default()
                    .fg(AyuTheme::ORANGE)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![
                Span::styled("• Molecule: ", Style::default().fg(AyuTheme::FG_SECONDARY)),
                Span::styled(
                    &self.config.moleculename,
                    Style::default().fg(AyuTheme::YELLOW),
                ),
            ]),
            Line::from(vec![
                Span::styled("• Sequence: ", Style::default().fg(AyuTheme::FG_SECONDARY)),
                Span::styled(&self.config.sequence, Style::default().fg(AyuTheme::GREEN)),
            ]),
            Line::from(vec![
                Span::styled("• Samples: ", Style::default().fg(AyuTheme::FG_SECONDARY)),
                Span::styled(
                    self.config.n_samples.to_string(),
                    Style::default().fg(AyuTheme::CYAN),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "• Top N for MD: ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    self.config.md_simulation.to_string(),
                    Style::default().fg(AyuTheme::PURPLE),
                ),
            ]),
            Line::from(vec![
                Span::styled("• Method: ", Style::default().fg(AyuTheme::FG_SECONDARY)),
                Span::styled(&self.config.method, Style::default().fg(AyuTheme::BLUE)),
            ]),
            Line::from(vec![
                Span::styled(
                    "• Temperature: ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    format!("{} K", self.config.temp),
                    Style::default().fg(AyuTheme::RED),
                ),
            ]),
            Line::from(vec![
                Span::styled("• Device: ", Style::default().fg(AyuTheme::FG_SECONDARY)),
                Span::styled(&self.config.device, Style::default().fg(AyuTheme::ORANGE)),
            ]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Parameter Info:",
                Style::default()
                    .fg(AyuTheme::GREEN)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![Span::styled(
                "• MD Simulations: Top N samples for",
                Style::default().fg(AyuTheme::FG_SECONDARY),
            )]),
            Line::from(vec![Span::styled(
                "  both MD simulation and final output",
                Style::default().fg(AyuTheme::FG_SECONDARY),
            )]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Available Methods:",
                Style::default()
                    .fg(AyuTheme::BLUE)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![Span::styled(
                "• fold - Protein folding simulation",
                Style::default().fg(AyuTheme::FG_SECONDARY),
            )]),
            Line::from(vec![Span::styled(
                "• search - Conformational search",
                Style::default().fg(AyuTheme::FG_SECONDARY),
            )]),
            Line::from(vec![Span::styled(
                "• explore - Exploration sampling",
                Style::default().fg(AyuTheme::FG_SECONDARY),
            )]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Available Devices:",
                Style::default()
                    .fg(AyuTheme::PURPLE)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![Span::styled(
                "• CPU - CPU computation",
                Style::default().fg(AyuTheme::FG_SECONDARY),
            )]),
            Line::from(vec![Span::styled(
                "• CUDA - NVIDIA GPU",
                Style::default().fg(AyuTheme::FG_SECONDARY),
            )]),
            Line::from(vec![Span::styled(
                "• HIP - AMD GPU",
                Style::default().fg(AyuTheme::FG_SECONDARY),
            )]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Press ENTER to edit selected field",
                Style::default().fg(AyuTheme::CYAN),
            )]),
            Line::from(vec![Span::styled(
                "Press F5 to save configuration",
                Style::default().fg(AyuTheme::YELLOW),
            )]),
        ];

        let details_paragraph = Paragraph::new(details)
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .title("Help & Info")
                    .border_style(Style::default().fg(AyuTheme::BORDER))
                    .title_style(
                        Style::default()
                            .fg(AyuTheme::PURPLE)
                            .add_modifier(Modifier::BOLD),
                    ),
            )
            .wrap(ratatui::widgets::Wrap { trim: true });
        f.render_widget(details_paragraph, chunks[1]);
    }

    fn render_configuration_review(&mut self, f: &mut Frame, area: ratatui::layout::Rect) {
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(5), // Header with selected peptide info
                Constraint::Min(0),    // Configuration fields
                Constraint::Length(4), // Action buttons
            ])
            .split(area);

        // Header with selected peptide info
        let header_text = vec![
            Line::from(vec![
                Span::styled(
                    "Selected Peptide: ",
                    Style::default()
                        .fg(AyuTheme::ORANGE)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    &self.config.moleculename,
                    Style::default()
                        .fg(AyuTheme::YELLOW)
                        .add_modifier(Modifier::BOLD),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "Sequence: ",
                    Style::default()
                        .fg(AyuTheme::GREEN)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(&self.config.sequence, Style::default().fg(AyuTheme::CYAN)),
            ]),
            Line::from(vec![
                Span::styled(
                    "Review and edit configuration below, then press ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    "'r'",
                    Style::default()
                        .fg(AyuTheme::RED)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    " to run simulation:",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
            ]),
        ];

        let header = Paragraph::new(header_text)
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .title("Configuration Review")
                    .border_style(Style::default().fg(AyuTheme::BORDER))
                    .title_style(
                        Style::default()
                            .fg(AyuTheme::CYAN)
                            .add_modifier(Modifier::BOLD),
                    ),
            )
            .style(Style::default().fg(AyuTheme::FG_PRIMARY));
        f.render_widget(header, chunks[0]);

        // Configuration fields (reuse the configuration rendering logic)
        let config_chunks = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Percentage(60), Constraint::Percentage(40)])
            .split(chunks[1]);

        // Configuration fields list
        let config_items: Vec<ListItem> = self
            .config_fields
            .iter()
            .map(|field| {
                let (name, value) = match field {
                    ConfigField::Moleculename => {
                        ("Molecule Name", self.config.moleculename.clone())
                    }
                    ConfigField::Folder => ("Output Folder", self.config.folder.clone()),
                    ConfigField::Sequence => (
                        "Sequence",
                        if self.config.sequence.len() > 25 {
                            format!("{}...", &self.config.sequence[..22])
                        } else {
                            self.config.sequence.clone()
                        },
                    ),
                    ConfigField::Interface => ("Interface ⚡", self.config.interface.clone()),
                    ConfigField::NSamples => {
                        ("Number of Samples", self.config.n_samples.to_string())
                    }
                    ConfigField::MdSimulation => (
                        "MD Simulations (Top N)",
                        self.config.md_simulation.to_string(),
                    ),
                    ConfigField::Temp => ("Temperature (K)", self.config.temp.to_string()),
                    ConfigField::Forcefield => (
                        "Force Fields ⚡",
                        if self.config.forcefield.len() <= 2 {
                            self.config.forcefield.join(", ")
                        } else {
                            format!(
                                "{} fields: {}, ...",
                                self.config.forcefield.len(),
                                &self.config.forcefield[0]
                            )
                        },
                    ),
                    ConfigField::Device => ("Device ⚡", self.config.device.clone()),
                    ConfigField::Solvent => ("Solvent Count", self.config.solvent.to_string()),
                    ConfigField::Steps => ("Equilibration Steps", self.config.steps.to_string()),
                    ConfigField::Gamma => ("Friction Coefficient", self.config.gamma.to_string()),
                    ConfigField::Dt => ("Time Step (ps)", self.config.dt.to_string()),
                    ConfigField::MdSteps => ("MD Steps", self.config.md_steps.to_string()),
                    ConfigField::Method => ("Sampling Method ⚡", self.config.method.clone()),
                };

                let content = if Some(field) == self.editing_field.as_ref() {
                    if self.enum_selection_mode {
                        format!("{:<25} │ ▼ {} ◄ SELECTING", name, self.edit_buffer)
                    } else {
                        format!("{:<25} │   {}_", name, self.edit_buffer)
                    }
                } else {
                    let has_enum_options = matches!(
                        field,
                        ConfigField::Interface
                            | ConfigField::Device
                            | ConfigField::Method
                            | ConfigField::Forcefield
                    );
                    if has_enum_options {
                        format!("{:<25} │ ▼ {}", name, value)
                    } else {
                        format!("{:<25} │   {}", name, value)
                    }
                };

                ListItem::new(content)
            })
            .collect();

        let config_list = List::new(config_items)
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .title("Configuration Parameters")
                    .border_style(Style::default().fg(AyuTheme::BORDER))
                    .title_style(
                        Style::default()
                            .fg(AyuTheme::ORANGE)
                            .add_modifier(Modifier::BOLD),
                    ),
            )
            .style(Style::default().fg(AyuTheme::FG_PRIMARY))
            .highlight_style(
                Style::default()
                    .bg(AyuTheme::SELECTION)
                    .fg(AyuTheme::BLUE)
                    .add_modifier(Modifier::BOLD),
            )
            .highlight_symbol("▶ ");

        f.render_stateful_widget(&config_list, config_chunks[0], &mut self.config_list_state);

        // Render enum selection popup if in enum mode
        if self.enum_selection_mode && !self.enum_options.is_empty() {
            self.render_enum_popup(f, config_chunks[0]);
        }

        // Configuration summary and validation
        let summary_text = vec![
            Line::from(vec![Span::styled(
                "Simulation Summary",
                Style::default()
                    .add_modifier(Modifier::BOLD)
                    .fg(AyuTheme::GREEN),
            )]),
            Line::from(""),
            Line::from(vec![
                Span::styled(
                    "• Will generate ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    self.config.n_samples.to_string(),
                    Style::default().fg(AyuTheme::CYAN),
                ),
                Span::styled(" samples", Style::default().fg(AyuTheme::FG_SECONDARY)),
            ]),
            Line::from(vec![
                Span::styled("• Using ", Style::default().fg(AyuTheme::FG_SECONDARY)),
                Span::styled(&self.config.method, Style::default().fg(AyuTheme::BLUE)),
                Span::styled(" method", Style::default().fg(AyuTheme::FG_SECONDARY)),
            ]),
            Line::from(vec![
                Span::styled("• Top ", Style::default().fg(AyuTheme::FG_SECONDARY)),
                Span::styled(
                    self.config.md_simulation.to_string(),
                    Style::default().fg(AyuTheme::PURPLE),
                ),
                Span::styled(
                    " samples for MD",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "• Temperature: ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    format!("{} K", self.config.temp),
                    Style::default().fg(AyuTheme::RED),
                ),
            ]),
            Line::from(vec![
                Span::styled("• Device: ", Style::default().fg(AyuTheme::FG_SECONDARY)),
                Span::styled(&self.config.device, Style::default().fg(AyuTheme::ORANGE)),
            ]),
            Line::from(vec![
                Span::styled(
                    "• Force fields: ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    self.config.forcefield.join(", "),
                    Style::default().fg(AyuTheme::YELLOW),
                ),
            ]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Estimated Runtime:",
                Style::default()
                    .fg(AyuTheme::GREEN)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![Span::styled(
                self.estimate_runtime(),
                Style::default().fg(AyuTheme::CYAN),
            )]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Output Location:",
                Style::default()
                    .fg(AyuTheme::PURPLE)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![Span::styled(
                format!("{}/{}/", self.config.folder, self.config.moleculename),
                Style::default().fg(AyuTheme::YELLOW),
            )]),
        ];

        let summary = Paragraph::new(summary_text)
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .title("Summary")
                    .border_style(Style::default().fg(AyuTheme::BORDER))
                    .title_style(
                        Style::default()
                            .fg(AyuTheme::GREEN)
                            .add_modifier(Modifier::BOLD),
                    ),
            )
            .wrap(ratatui::widgets::Wrap { trim: true });
        f.render_widget(summary, config_chunks[1]);

        // Action buttons
        let action_text = vec![
            Line::from(""),
            Line::from(vec![
                Span::styled("Press ", Style::default().fg(AyuTheme::FG_SECONDARY)),
                Span::styled(
                    "'r'",
                    Style::default()
                        .fg(AyuTheme::GREEN)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    " to RUN SIMULATION with current settings",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
            ]),
            Line::from(vec![
                Span::styled("Press ", Style::default().fg(AyuTheme::FG_SECONDARY)),
                Span::styled(
                    "ESC",
                    Style::default()
                        .fg(AyuTheme::RED)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    " to go back to peptide selection",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
            ]),
        ];

        let actions = Paragraph::new(action_text)
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .title("Actions")
                    .border_style(Style::default().fg(AyuTheme::BORDER))
                    .title_style(
                        Style::default()
                            .fg(AyuTheme::RED)
                            .add_modifier(Modifier::BOLD),
                    ),
            )
            .alignment(Alignment::Center);
        f.render_widget(actions, chunks[2]);
    }

    fn estimate_runtime(&self) -> String {
        // Simple runtime estimation based on samples and method
        let base_time = match self.config.method.as_str() {
            "fold" => 2.0,
            "search" => 1.5,
            "explore" => 3.0,
            _ => 2.0,
        };

        let sample_factor = (self.config.n_samples as f64 / 1000.0).max(0.1);
        let md_factor = (self.config.md_simulation as f64 / 10.0).max(0.1);

        let estimated_minutes = base_time * sample_factor * md_factor;

        if estimated_minutes < 1.0 {
            format!("~{:.0} seconds", estimated_minutes * 60.0)
        } else if estimated_minutes < 60.0 {
            format!("~{:.1} minutes", estimated_minutes)
        } else {
            format!("~{:.1} hours", estimated_minutes / 60.0)
        }
    }

    fn render_enum_popup(&mut self, f: &mut Frame, area: ratatui::layout::Rect) {
        // Ensure we have options to display
        if self.enum_options.is_empty() {
            return;
        }

        // Calculate popup size based on options
        let max_option_length = self
            .enum_options
            .iter()
            .map(|opt| opt.len())
            .max()
            .unwrap_or(20);

        let popup_height = (self.enum_options.len() + 5).min(15) as u16; // +5 for borders and instructions
        let popup_width = (max_option_length + 12).max(35) as u16; // +12 for padding and indicators

        // Center the popup with better positioning
        let popup_area = centered_rect(
            ((popup_width * 100) / area.width).min(85),
            ((popup_height * 100) / area.height).min(75),
            area,
        );

        // Split popup area for list and instructions
        let popup_chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Min(0),    // List area
                Constraint::Length(3), // Instructions area
            ])
            .split(popup_area);

        // Create enum option list items with enhanced formatting
        let enum_items: Vec<ListItem> = self
            .enum_options
            .iter()
            .enumerate()
            .map(|(i, option)| {
                let is_selected = i == self.enum_selected_index;
                let content = if is_selected {
                    format!("► {}", option)
                } else {
                    format!("  {}", option)
                };
                ListItem::new(content).style(if is_selected {
                    Style::default()
                        .fg(AyuTheme::CYAN)
                        .add_modifier(Modifier::BOLD)
                } else {
                    Style::default().fg(AyuTheme::FG_PRIMARY)
                })
            })
            .collect();

        // Get field name for title
        let field_name = if let Some(field) = &self.editing_field {
            match field {
                ConfigField::Interface => "Interface Options",
                ConfigField::Device => "Compute Device",
                ConfigField::Method => "Sampling Method",
                ConfigField::Forcefield => "Force Field",
                _ => "Available Options",
            }
        } else {
            "Available Options"
        };

        let enum_list = List::new(enum_items)
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .title(format!("▼ {} (Use ↑/↓ to navigate)", field_name))
                    .border_style(Style::default().fg(AyuTheme::YELLOW))
                    .title_style(
                        Style::default()
                            .fg(AyuTheme::YELLOW)
                            .add_modifier(Modifier::BOLD),
                    ),
            )
            .style(Style::default().bg(AyuTheme::BG_SECONDARY))
            .highlight_style(
                Style::default()
                    .fg(AyuTheme::GREEN)
                    .bg(AyuTheme::SELECTION)
                    .add_modifier(Modifier::BOLD),
            )
            .highlight_symbol("");

        // Clear the area and render the popup
        f.render_widget(Clear, popup_area);
        f.render_stateful_widget(enum_list, popup_chunks[0], &mut self.enum_list_state);

        // Add instructions in a separate bordered area
        let instructions = Paragraph::new(vec![Line::from(vec![
            Span::styled(
                "↑/↓",
                Style::default()
                    .fg(AyuTheme::YELLOW)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::styled(" Navigate  ", Style::default().fg(AyuTheme::FG_SECONDARY)),
            Span::styled(
                "ENTER",
                Style::default()
                    .fg(AyuTheme::GREEN)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::styled(
                " Select Option",
                Style::default().fg(AyuTheme::FG_SECONDARY),
            ),
        ])])
        .block(
            Block::default()
                .borders(Borders::ALL)
                .border_style(Style::default().fg(AyuTheme::BORDER)),
        )
        .alignment(Alignment::Center)
        .style(Style::default().bg(AyuTheme::BG_SECONDARY));

        f.render_widget(instructions, popup_chunks[1]);
    }

    fn render_help(&self, f: &mut Frame) {
        let help_text = vec![
            Line::from(vec![Span::styled(
                "DNCS Enhanced Control Panel - Help",
                Style::default()
                    .add_modifier(Modifier::BOLD)
                    .fg(AyuTheme::CYAN),
            )]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Global Keys:",
                Style::default()
                    .fg(AyuTheme::ORANGE)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![
                Span::styled(
                    "  TAB         - ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    "Switch between Peptide Selection and Configuration",
                    Style::default().fg(AyuTheme::FG_PRIMARY),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "  F5          - ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    "Save current configuration to dncs.toml",
                    Style::default().fg(AyuTheme::FG_PRIMARY),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "  F9          - ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    "Run simulation with current settings",
                    Style::default().fg(AyuTheme::FG_PRIMARY),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "  h or F1     - ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    "Toggle this help screen",
                    Style::default().fg(AyuTheme::FG_PRIMARY),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "  q or ESC    - ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    "Quit application",
                    Style::default().fg(AyuTheme::FG_PRIMARY),
                ),
            ]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Peptide Selection Mode:",
                Style::default()
                    .fg(AyuTheme::GREEN)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![
                Span::styled(
                    "  ↑/↓ or j/k  - ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    "Navigate peptide list",
                    Style::default().fg(AyuTheme::FG_PRIMARY),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "  ENTER       - ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    "Select peptide → Review configuration",
                    Style::default().fg(AyuTheme::YELLOW),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "  /           - ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled("Start search", Style::default().fg(AyuTheme::FG_PRIMARY)),
            ]),
            Line::from(vec![
                Span::styled(
                    "  s           - ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    "Cycle sort modes (Length/Name/Type)",
                    Style::default().fg(AyuTheme::FG_PRIMARY),
                ),
            ]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Configuration Review Mode:",
                Style::default()
                    .fg(AyuTheme::BLUE)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![
                Span::styled(
                    "  r           - ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    "RUN SIMULATION (just run)",
                    Style::default()
                        .fg(AyuTheme::GREEN)
                        .add_modifier(Modifier::BOLD),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "  ENTER       - ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    "Edit selected field",
                    Style::default().fg(AyuTheme::FG_PRIMARY),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "  ESC         - ",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
                Span::styled(
                    "Back to peptide selection",
                    Style::default().fg(AyuTheme::FG_PRIMARY),
                ),
            ]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Smart Field Editing:",
                Style::default()
                    .fg(AyuTheme::ORANGE)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![Span::styled(
                "  • Fields with ⚡ show dropdown menus",
                Style::default().fg(AyuTheme::CYAN),
            )]),
            Line::from(vec![Span::styled(
                "  • Use ↑/↓ to select from predefined options",
                Style::default().fg(AyuTheme::CYAN),
            )]),
            Line::from(vec![Span::styled(
                "  • Numeric fields allow direct typing",
                Style::default().fg(AyuTheme::YELLOW),
            )]),
            Line::from(vec![Span::styled(
                "  • 🔽 indicates active dropdown selection",
                Style::default().fg(AyuTheme::PURPLE),
            )]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Status Indicators:",
                Style::default()
                    .fg(AyuTheme::GREEN)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![
                Span::styled(
                    "  ✓ ",
                    Style::default()
                        .fg(AyuTheme::GREEN)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    "Results available for this molecule",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
            ]),
            Line::from(vec![
                Span::styled("  ○ ", Style::default().fg(AyuTheme::FG_MUTED)),
                Span::styled(
                    "Molecule not yet processed",
                    Style::default().fg(AyuTheme::FG_SECONDARY),
                ),
            ]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Sorting Options:",
                Style::default()
                    .fg(AyuTheme::PURPLE)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![Span::styled(
                "  Press 's' to cycle through:",
                Style::default().fg(AyuTheme::FG_SECONDARY),
            )]),
            Line::from(vec![Span::styled(
                "  • Length ↑ - Shortest to longest",
                Style::default().fg(AyuTheme::CYAN),
            )]),
            Line::from(vec![Span::styled(
                "  • Length ↓ - Longest to shortest",
                Style::default().fg(AyuTheme::CYAN),
            )]),
            Line::from(vec![Span::styled(
                "  • Name ↑ - Alphabetical A-Z",
                Style::default().fg(AyuTheme::YELLOW),
            )]),
            Line::from(vec![Span::styled(
                "  • Name ↓ - Alphabetical Z-A",
                Style::default().fg(AyuTheme::YELLOW),
            )]),
            Line::from(vec![Span::styled(
                "  • Type ↑/↓ - By peptide type",
                Style::default().fg(AyuTheme::GREEN),
            )]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Directory Management:",
                Style::default()
                    .fg(AyuTheme::RED)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![Span::styled(
                "  When existing sample directories are found:",
                Style::default().fg(AyuTheme::FG_SECONDARY),
            )]),
            Line::from(vec![Span::styled(
                "  • [D] Delete and continue",
                Style::default().fg(AyuTheme::RED),
            )]),
            Line::from(vec![Span::styled(
                "  • [B] Backup with timestamp",
                Style::default().fg(AyuTheme::YELLOW),
            )]),
            Line::from(vec![Span::styled(
                "  • [R] Rename with timestamp",
                Style::default().fg(AyuTheme::BLUE),
            )]),
            Line::from(vec![Span::styled(
                "  • [C] Cancel operation",
                Style::default().fg(AyuTheme::FG_MUTED),
            )]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "Press any key to return to the main interface...",
                Style::default()
                    .fg(AyuTheme::ORANGE)
                    .add_modifier(Modifier::ITALIC),
            )]),
        ];

        let help_paragraph = Paragraph::new(help_text)
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .title("Help")
                    .border_style(Style::default().fg(AyuTheme::CYAN))
                    .title_style(
                        Style::default()
                            .fg(AyuTheme::CYAN)
                            .add_modifier(Modifier::BOLD),
                    ),
            )
            .style(
                Style::default()
                    .bg(AyuTheme::BG_SECONDARY)
                    .fg(AyuTheme::FG_PRIMARY),
            )
            .alignment(Alignment::Left)
            .wrap(ratatui::widgets::Wrap { trim: true });

        let area = centered_rect(80, 90, f.area());
        f.render_widget(Clear, area);
        f.render_widget(help_paragraph, area);
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
